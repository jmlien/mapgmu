#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/time.h> //fps calculations
#include <pcl/common/transforms.h>
#include <pcl/io/hdl_grabber.h>
#include <pcl/visualization/point_cloud_color_handlers.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/image_viewer.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/boost.h>
#include <pcl/visualization/mouse_event.h>
#include <vector>
#include <string>
#include <Eigen/Eigen>
#include <typeinfo>
#include <pcap/pcap.h>
#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/asio.hpp>
#include <boost/version.hpp>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/math/special_functions.hpp>
#include <pcl/features/don.h>
#include <pcl/search/organized.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/filters/extract_indices.h>
#include <stdio.h>
#include "MyMeshBuilder.h"
#include "DonEstimator.hpp"
#include "mapgmuutils.hpp"
#include "CubeAligner.hpp"

using namespace std;
using namespace pcl;
using namespace pcl::console;
using namespace pcl::visualization;
using namespace boost::filesystem;

typedef pcl::PointNormal PointNT;
typedef pcl::PointCloud<PointNT> NCloud;
typedef pcl::PointCloud<PointNT>::Ptr NCloudPtr;
typedef CubeAligner::CCluster Cluster;

double startTime = -1;
char timeFile[100];
int _min = -1;
int _max = -1;
int _skip = -1;
string _kml_file = "";
mapgmu::PathFinder *_pfp = 0;

template<typename PointType>
class SimplePCDViewer
{
  public:
    typedef PointCloud<PointType> Cloud;
    typedef typename Cloud::Ptr CloudPtr;
    typedef typename Cloud::ConstPtr CloudConstPtr;

    SimplePCDViewer (std::string file,
                     PointCloudColorHandler<PointType> &handler,
		     bool isDir, bool hasNormals, bool meshed, bool transformFile
		    )
      : fileName_ (file)
      , handler_ (handler)
      , isDir_(isDir)
      , hasNorms_(hasNormals)
      , meshed_(meshed)
      , transfile_(transformFile)
      , quitFlag(false)
      , color_handler_int("intensity")
      , color_handler_curv("curvature")
      , isIntOn (false)
      , isSpheresOn (false)
      , cloud_ (new Cloud())
      , cpcloud_ (new Cloud())
      , ncloud_ (new NCloud())
      , cloud_viewer_ (new PCLVisualizer ("PCL HDL Cloud"))
    {

    }

    void LoadTransformFile()
    {
      FILE *fp;
      Eigen::Matrix4d transMat = Eigen::Matrix4d::Identity();
      Eigen::Matrix4d tempMat;
      Cloud tmpCloud;
      typename Cloud::iterator cloudItr;
      CloudPtr transedCloud;
      std::string tmpStr;
      char myFName[50];
      double tmpFloat;
      int i,j,k=0;
      ExtractIndices<PointXYZI> iextract(false);
      mapgmu::PathFinder pf;
      PointXYZ p1,p2;
      p1.x=0;p1.y=0;p1.z=0;


      fp = fopen(fileName_.c_str(),"r");
      cout << "Reading Transform File " << fileName_ << "..." << endl;

      while (fscanf(fp,"%s",myFName) != EOF)
      {
      	if (strlen(myFName) == 0){break;}
      	for (i=0;i<4;i++){for(j=0;j<4;j++){
      	  fscanf(fp,"%le",&tmpFloat);
      	  tempMat(i,j) = tmpFloat;
      	}}
      	k++;
      	transMat = transMat*tempMat;
      	p2 = p1;
      	p1 = pf.extractOriginD(transMat);
      	tmpStr = std::string(myFName);
      	sprintf(myFName,"Line%d",k);
      	//add a line to show the path of the clouds
      	cloud_viewer_->addLine(p1,p2,
      			       ((k*33)%255)/255.0,
      			       ((k*66)%255)/255.0,
      			       ((k*99)%255)/255.0,std::string(myFName),0);

      	//only add the clouds that are not missed by skipping and min/max
      	if ((_skip < 0 || (k%_skip)==0) &&
      	  (_min < 0 || k >= _min) &&
      	  (_max <  0 || k <= _max) )
      	{
      	  cout << "Loading pcd file " << tmpStr << "..." << endl;
      	  transedCloud = CloudPtr(new Cloud);
      	  pcl::io::loadPCDFile(tmpStr.c_str(),tmpCloud);
      	  pcl::transformPointCloud (tmpCloud,*transedCloud,transMat);

      	  clouds.push_back(transedCloud);
      	  cloudNames.push_back(tmpStr);
      	  handler_.setInputCloud(transedCloud);
      	    cloud_viewer_->addPointCloud(transedCloud,handler_,tmpStr);
      	}
      }

      fclose(fp);

    }//end of LoadTransformFile

    void LoadKmlFile(string kmlFileName)
    {
      _pfp = new mapgmu::PathFinder(kmlFileName,"my path");
      int segCount = _pfp->getWaypointCount();
      int i;
      char myLine[10];

      cout << "Loaded " << segCount << " waypoints" << endl;

      for (i=1;i<segCount;i++)
      {
      	//add a line segment for each waypoint after the first
      	sprintf(myLine,"l%d",i);
      	cloud_viewer_->addLine(_pfp->getWaypoint(i),
			                         _pfp->getWaypoint(i-1),
			                         ((i*33)%255)/255.0,
			                         ((i*66)%255)/255.0,
			                         ((i*99)%255)/255.0,
			                         std::string(myLine),0);
      }//end for i
    }//LoadKmlFile


    void keyboard_callback(const pcl::visualization::KeyboardEvent &event)
    {
      typename std::vector<CloudPtr>::iterator cloudItr;
      typename Cloud::iterator pointItr;
      int i=0,rIdx;
      double sRad;
      char sName[10];
      float f;
      Eigen::Matrix4f timeTransMat;
      CloudPtr transCloud(new Cloud);

      //only trigger on key_down
      if (!event.keyUp())
      {
	       return;
      }

      //quit
      if (event.getKeyCode() == 'q')
      {
	       quitFlag = true;
	       return;
      }

      //toggle small spheres, size of DoNFilters getMinRadiusForPoint
      if (event.getKeyCode() == 's')
      {
      	if (!isSpheresOn)
      	{
      	  for (cloudItr = clouds.begin();cloudItr!=clouds.end();cloudItr++)
      	  {
      	    for (i=0;i<100;i++)
      	    {
      	      rIdx = rand() % (*cloudItr)->size();
      	      sRad = donEst.getMinRadiusForPoint((*cloudItr)->at(rIdx));
      	      sprintf(sName,"%d_%d",i,rIdx);
      	      cloud_viewer_->addSphere((*cloudItr)->at(rIdx),sRad,sName);
      	    }
      	  }
      	}
      	else
      	{
      	  cloud_viewer_->removeAllShapes();
      	}
      	isSpheresOn = !isSpheresOn;
      	return;
      }

      //Color each cloud with different color,
      if (event.getKeyCode() == 'L')
      {
  	    cloud_viewer_->removeAllPointClouds();
      	for (cloudItr = clouds.begin();cloudItr!=clouds.end();cloudItr++)
      	{
      	  PointCloudColorHandlerCustom<PointXYZI> cHandNew(230-(i*33)%230,220-(i*66)%220,240-(i*99)%240);
      	  cHandNew.setInputCloud(*cloudItr);
      	  if (!cloud_viewer_->updatePointCloud(*cloudItr,cHandNew,cloudNames.at(i)))
      	    cloud_viewer_->addPointCloud(*cloudItr,cHandNew,cloudNames.at(i));
      	  i++;
      	}
      	return;
      }

      //toggle large spheres, size of DoNFilter's getMaxRadiusForPoint
      if (event.getKeyCode() == 'S')
      {
        	if (!isSpheresOn)
        	{
        	  for (cloudItr = clouds.begin();cloudItr!=clouds.end();cloudItr++)
        	  {
        	    for (i=0;i<100;i++)
        	    {
        	      rIdx = rand() % (*cloudItr)->size();
        	      sRad = donEst.getMaxRadiusForPoint((*cloudItr)->at(rIdx));
        	      sprintf(sName,"%d_%d",i,rIdx);
        	      cloud_viewer_->addSphere((*cloudItr)->at(rIdx),sRad,sName);
        	    }
        	  }
        	}
        	else
        	{
        	  cloud_viewer_->removeAllShapes();
        	}
        	isSpheresOn = !isSpheresOn;
        	return;
      }

      //toggle coloring by intensity
      if (event.getKeyCode() == 'i')
      {
        	if (isIntOn)
        	{
        	  if (!hasNorms_)
        	  {
        	    for (cloudItr = clouds.begin();cloudItr!=clouds.end();cloudItr++)
        	    {
        	      handler_.setInputCloud(*cloudItr);
        	      if (!cloud_viewer_->updatePointCloud(*cloudItr,handler_,cloudNames.at(i)))
        		cloud_viewer_->addPointCloud(*cloudItr,handler_,cloudNames.at(i));
        	      i++;
        	    }
        	  }
        	  else
        	  {
        	    cloud_viewer_->removeAllPointClouds();
        	    cloud_viewer_->addPointCloudNormals<pcl::PointNormal>(ncloud_,1,0.1,"DoN Cloud",0);
        	  }
        	}
        	else
        	{
        	  if (!hasNorms_)
        	  {
        	    for (cloudItr = clouds.begin();cloudItr!=clouds.end();cloudItr++)
        	    {
        	      color_handler_int.setInputCloud(*cloudItr);
        	      if (!cloud_viewer_->updatePointCloud(*cloudItr,color_handler_int,cloudNames.at(i)))
        		cloud_viewer_->addPointCloud(*cloudItr,color_handler_int,cloudNames.at(i));
        	      i++;
        	    }
        	  }
        	  else
        	  {
        	    color_handler_curv.setInputCloud(ncloud_);
        	    cloud_viewer_->removeAllPointClouds();
        	    cloud_viewer_->addPointCloud(ncloud_,color_handler_curv,"DoN Cloud");

        	  }
        	}
        	isIntOn = !isIntOn;
        	return;
      }

      //toggle polygons from clusters created by CubeAligner
      /*if (event.getKeyCode() == 'u' )
      {
      	CubeAligner cAlign(cloud_,0.5);

      	cAlign.identifyClusters();

      	std::vector<Cluster> *clusters = cAlign.getClusters();

      	std::vector<Cluster>::iterator clIt;
      	i=0;

      	for (clIt = clusters->begin(); clIt < clusters->end();
      	  clIt++)
      	{
      	  //if (i>3)break;
      	  sprintf(sName,"cluster_%d",++i);
      	  cloud_viewer_->removeShape(sName);
      	  cloud_viewer_->addPolygon<PointType>((*clIt).toPolygon(i),1.0,1.0,0.0,std::string(sName),0);

      	  sprintf(sName,"sphere_%d",i,rIdx);
      	  cloud_viewer_->removeShape(sName);
      	  cloud_viewer_->addSphere((*clIt).getCOG(),1.0,sName);
      	}

      	return;
      }*/

      //advance or reverse clouds if working with directory
      if (isDir_)
      {
      	//goes backwards one cloud
      	if (event.getKeyCode() == 'b')
      	{
      	  fileIter--;
      	  if (fileIter->extension() == ".ts")
      	  {
      	    fileIter-=3;
      	  }
      	  else
      	  {
      	    fileIter--;
      	  }
      	}
      	//
      	if (!cloud_viewer_->updateText("",200,0,26,255,255,255,"text"))
      	  cloud_viewer_->addText("",200,0,26,255,255,255,"text");
      	// See if we can get a cloud
      	if (hasNorms_)
      	{
      	  if (!pcl::io::loadPCDFile(fileIter->c_str(),*ncloud_))
      	  {
      	    if (isIntOn)
      	    {
      	      color_handler_curv.setInputCloud(ncloud_);
      	      cloud_viewer_->removeAllPointClouds();
      	      cloud_viewer_->addPointCloud(ncloud_,color_handler_curv,"DoN Cloud");
      	    }
      	    else
      	    {
      	      cloud_viewer_->removeAllPointClouds();
      	      cloud_viewer_->addPointCloudNormals<pcl::PointNormal>(ncloud_,1,0.1,"DoN Cloud",0);
      	    }
      	    //display the cloud name
      	    if (!cloud_viewer_->updateText(fileIter->c_str(),200,0,26,1.0,1.0,1.0,"text"))
      	      cloud_viewer_->addText(fileIter->c_str(),200,0,26,1.0,1.0,1.0,"text");
      	  }
      	}
      	else
      	{
      	  if (!pcl::io::loadPCDFile(fileIter->c_str(),*cloud_))
      	  {
      	    if (fileIter > files.begin())
      	    {
      	      CubeAligner cb(cpcloud_,0.1);
      	      cpcloud_->clear();
      	      *cpcloud_+=*cloud_;
      	      //display cubic filtered cloud with 'c'
      	      if (event.getKeyCode() == 'c')
      		cb.filterCloudCubic(cloud_);
      	    }
      	    else
      	    {
      	      cpcloud_->clear();
      	      *cpcloud_+=*cloud_;
      	    }
      	    fileIter++;

      	    //parse the time file if available
      	    if (fileIter->extension() == ".ts")
      	    {
      	      FILE *tsFile = fopen(fileIter->c_str(),"r");
      	      double currTime;
      	      fscanf(tsFile,"%lf",&currTime);
      	      fclose(tsFile);
      	      if (startTime == -1)
      	      {
      		startTime = currTime;
      	      }
      	      sprintf(timeFile,"%.4f sec _ %s",(currTime-startTime)*1.0e-6,fileIter->c_str());
      	      //if we have a path, move the cloud along it
      	      if (_pfp)
      	      {
      		timeTransMat = _pfp->getTransForTime((currTime-startTime)*1.0e-6);
      		transCloud->clear();
      		pcl::transformPointCloud (*cloud_,*transCloud,timeTransMat);
      		cloud_->clear();
      		*cloud_+=*transCloud;
      	      }
      	    }
      	    else
      	    {
      	      fileIter--;
      	      sprintf(timeFile,"%s",fileIter->c_str());
      	    }

      	    //display cloud with handler based on intensity setting
      	    if (isIntOn)
      	    {
      	      color_handler_int.setInputCloud(cloud_);
      	      if (!cloud_viewer_->updatePointCloud(cloud_,color_handler_int,"PCD Cloud"))
      		cloud_viewer_->addPointCloud(cloud_,color_handler_int,"PCD Cloud");
      	    }
      	    else
      	    {
      	      handler_.setInputCloud(cloud_);
      	      if (!cloud_viewer_->updatePointCloud(cloud_,handler_,"PCD Cloud"))
      		cloud_viewer_->addPointCloud(cloud_,handler_,"PCD Cloud");
      	    }
      	    if (meshed_)
      	    {
      	      std::vector<pcl::Vertices> meshVerts = donEst.buildMesh(cloud_);

      	      if (!cloud_viewer_->updatePolygonMesh<PointType>(cloud_,meshVerts,"Polygon"))
      		cloud_viewer_->addPolygonMesh<PointType>(cloud_,meshVerts,"Polygon",0);
      	    }

      	    //print the time text on the screen
      	    if (!cloud_viewer_->updateText(timeFile,200,0,26,1.0,1.0,1.0,"text"))
      	      cloud_viewer_->addText(timeFile,200,0,26,1.0,1.0,1.0,"text");
      	  }
      	}
      	fileIter++;
      }
    }

    void
    run ()
    {
      PointType p;

      donEst.setMeshFlag(true);

      cloud_viewer_->setBackgroundColor (0, 0, 0);

      //aerial view
      cloud_viewer_->setCameraPosition (3.15783,-25.694,73.0933, -1.43662,-17.8274,-1.5715, 0.146305,-0.982816,-0.112552, 0);
      //street view
      //cloud_viewer_->setCameraPosition (-43.2621,43.7286,12.0023,-0.598113,-13.1002,-12.6573,0.223411,-0.241993,0.944207,0);

      cloud_viewer_->setCameraClipDistances (0.0, 50.0);

      //build the file array if using directory
      if (isDir_)
      {
      	copy(directory_iterator(path(fileName_)),directory_iterator(),back_inserter(files));
      	std::sort(files.begin(),files.end());
      	fileIter = files.begin();
      	if (_min > 0)
      	{
      	  char fileNum[9];
      	  sprintf(fileNum,"%04d.pcd",_min);
      	  string fNumStr(fileNum);
      	  while (fileIter->string().find(fNumStr) == string::npos)
      	  {
      	    fileIter = files.erase(fileIter);
      	  }
      	  fileIter = files.begin();
      	}
      }
      // load the transform file and clouds
      else if (transfile_)
      {
	        LoadTransformFile();
      }
      //single cloud
      else
      {
      	if (hasNorms_){
      	  pcl::io::loadPCDFile(fileName_,*ncloud_);
      	  cloud_viewer_->addPointCloudNormals<pcl::PointNormal>(ncloud_,1,0.1,"DoN Cloud",0);
      	}
      	else
      	{
      	  pcl::io::loadPCDFile(fileName_,*cloud_);
      	  clouds.push_back(cloud_);
      	  cloudNames.push_back(fileName_);
      	  handler_.setInputCloud(cloud_);
      	    cloud_viewer_->addPointCloud(cloud_,handler_,fileName_);
      	  //add mesh if enabled
      	  if (meshed_)
      	  {
      	    std::vector<pcl::Vertices> meshVerts = donEst.buildMesh(cloud_);

      	    if (!cloud_viewer_->updatePolygonMesh<PointType>(cloud_,meshVerts,"Polygon"))
      	      cloud_viewer_->addPolygonMesh<PointType>(cloud_,meshVerts,"Polygon",0);
      	  }
      	}
      }

      cloud_viewer_->initCameraParameters ();

      //this crashes on OSX
      //cloud_viewer_->addCoordinateSystem (3.0);

      if (!_kml_file.empty())
      {
      	LoadKmlFile(_kml_file.c_str());
      	cout << "Loading KML File " << _kml_file << endl;
      }

      boost::function< void(const pcl::visualization::KeyboardEvent &)> cb = boost::bind(&SimplePCDViewer<PointType>::keyboard_callback,this,_1);
      cloud_viewer_->registerKeyboardCallback (cb);


      while (!quitFlag && (!isDir_ || fileIter!=files.end()))
      {
	       cloud_viewer_->spinOnce(100);
        //boost::this_thread::sleep (boost::posix_time::microseconds (100));
      }

    }//end of run


    boost::shared_ptr<PCLVisualizer> cloud_viewer_;
    std::string fileName_;
    PointCloudColorHandler<PointType> &handler_;
    bool isDir_;
    bool hasNorms_;
    bool meshed_;
    bool transfile_;
    bool quitFlag;
    CloudPtr cloud_;
    CloudPtr cpcloud_;
    NCloudPtr ncloud_;

    vector<path> files;
    vector<path>::iterator fileIter;
    DonEstimator<PointType> donEst;

    vector<CloudPtr> clouds;
    vector<std::string> cloudNames;

    PointCloudColorHandlerGenericField<PointXYZI> color_handler_int;
    PointCloudColorHandlerGenericField<PointNT> color_handler_curv;
    bool isIntOn;
    bool isSpheresOn;
};

void
usage (char ** argv)
{
  cout << "usage: " << argv[0] << endl
      << "\t-pcd <path-to-pcd-file>" << endl
      << "\t\t-OR-" << endl
      << "\t-dir <path-to-dir-containing-pcd-files>" << endl
      << "\t\t-OR-" << endl
      << "\t-transFile <path-to-transFile-output-from-cloud-reg>" << endl
      << "\t[-norm] (required to display clouds with normals)" << endl
      << "\t[-mesh] (creat mesh of each cloud, only for dir or pcd)" << endl
      << "\t[-min <minimum-file-index>] (dir or transFile only)" << endl
      << "\t[-max <maximum-file-index>] (transFile only)" << endl
      << "\t[-skip <use-nth-files>] (transFile only, eg 2 displays every other file)" << endl
      << "\t[-kml <path-to-kml-file>] (for dir, clouds will move along gps path)" << endl
      << endl;
  cout << argv[0] << " -h | --help : shows this help" << endl;
  return;
}

int
main (int argc, char ** argv)
{
  std::string pcdFile;

  if (
    find_switch (argc, argv, "-h") ||
    find_switch (argc, argv, "--help") ||
      !( find_switch (argc, argv, "-pcd") ||
	find_switch (argc, argv, "-dir") ||
	find_switch (argc, argv, "-transFile") ) ||
    (( find_switch (argc, argv, "-pcd") +
	find_switch (argc, argv, "-dir") +
	find_switch (argc, argv, "-transFile") ) > 1)
    )
  {
    usage (argv);
    return (0);
  }

  if (find_switch (argc, argv, "-pcd"))
  {
    parse_argument (argc, argv, "-pcd", pcdFile);
  }
  else if (find_switch (argc, argv, "-dir"))
  {
    if (find_switch (argc, argv, "-min") )
    {
      parse_argument (argc, argv, "-min", pcdFile);
      _min = atoi(pcdFile.c_str());
    }
    parse_argument (argc, argv, "-dir", pcdFile);
  }
  else if (find_switch (argc, argv, "-transFile"))
  {
    if (find_switch (argc, argv, "-min") )
    {
      parse_argument (argc, argv, "-min", pcdFile);
      _min = atoi(pcdFile.c_str());
    }
    if (find_switch (argc, argv, "-max") )
    {
      parse_argument (argc, argv, "-max", pcdFile);
      _max = atoi(pcdFile.c_str());
    }
    if (find_switch (argc, argv, "-skip") )
    {
      parse_argument (argc, argv, "-skip", pcdFile);
      _skip = atoi(pcdFile.c_str());
    }
    parse_argument (argc, argv, "-transFile", pcdFile);
    cout << "Parsed arg " << pcdFile << "...\n";
  }
  if (find_switch (argc, argv, "-kml"))
  {
    parse_argument (argc, argv, "-kml", _kml_file);
  }

  PointCloudColorHandlerCustom<PointXYZI> color_handler(255.0,255.0,255.0);

  SimplePCDViewer<PointXYZI> v (pcdFile,
				color_handler,
				find_switch (argc, argv, "-dir"),
				find_switch (argc, argv, "-norm"),
				find_switch (argc, argv, "-mesh"),
				find_switch (argc, argv, "-transFile")
			       );

  v.run ();
  return (0);
}
