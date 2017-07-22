#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/hdl_grabber.h>
#include <pcl/console/parse.h>
#include <string>
#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>
#include "MyHdlGrabber.h"
#include "DoNFilter.hpp"
#include "DonEstimator.hpp"
#include "CubeAligner.hpp"

using namespace std;
using namespace pcl;
using namespace pcl::console;
using namespace Eigen;


typedef PointCloud<PointXYZI> Cloud;
typedef PointCloud<PointXYZI>::Ptr CloudPtr;
typedef PointCloud<PointNormal> NCloud;
typedef PointCloud<PointNormal>::Ptr NCloudPtr;

void
usage (char ** argv)
{
  cout << "usage: " << argv[0] << endl
      << "\t-cal <path-to-calibration-file>" << endl
      << "\t-pcap <path-to-pcap-file>" << endl
      << "\t[-s <N>] will start at the Nth cloud" << endl
      << "\t[-don] output Density of Normal filtered clouds files in an inclosed dir */don" << endl
      << "\t[-intFilt <min-int-value>] save intensity filtered clouds to an inclosed dir */int" << endl
      << "\t[-cube] output a cubic polygon pcd file to */cube" << endl
      << "\t[-noCloud] stop the app from saving the unfiltered clouds, requires don,cube, or intFilt" << endl
      << "\t[-zFix] rotate the coordinate system to a vertical Z-axis (assuming velodyne MapGMU configuration)"
      << endl << endl
       << "Clouds are saved to <path-to-pcap-file(without '.pcap')>/####.pcd"
       << endl
       << "Timestamps are saved to <path-to-pcap-file(without '.pcap')>/####.ts" << endl << endl;
       
  cout << argv[0] << " -h | --help : shows this help" << endl;
  return;
}

int 
main (int argc, char ** argv)
{
  CloudPtr cloud_;
  CloudPtr intcloud_(new Cloud);
  Cloud aligned;
  NCloudPtr donCloud(new NCloud);
  std::string hdlCalibration, pcapFile,donDir,dirName,skipIdxStr,intDir,intStr,cubeDir;
  int cloudIdx =0;
  char fName[50];
  char fName2[50];
  FILE *tStampFile;
  DoNFilter<PointXYZI> donFilter;
  DonEstimator<PointXYZI> donEst;
  bool don,noCloud,intFilt,shouldRot,cube;
  int skipIdx = 0;
  double intScale;
  Eigen::Matrix4f transMat = Eigen::Matrix4f::Identity();
  Matrix3f rotMat;
  int i,j;
  Vector3f crossVec = Vector3f::Zero();
  double theta = -(3.141592/2.0);
  crossVec(0) = 1.0;
  rotMat = AngleAxisf(theta,crossVec);
  for(i=0;i<3;i++){for(j=0;j<3;j++){
    transMat(i,j) = rotMat(i,j);
  }}

  if (find_switch (argc, argv, "-h") || 
      find_switch (argc, argv, "--help") ||
      !(find_switch (argc, argv, "-cal") && find_switch (argc, argv, "-pcap")) ||
      (find_switch (argc, argv, "-noCloud") && !((find_switch (argc, argv, "-don") || find_switch(argc,argv,"-cube") ) ||
	find_switch (argc, argv, "-intFilt") )
      ))
  {
    usage (argv);
    return (0);
  }
  
  if (find_switch (argc, argv, "-s"))
  {
    parse_argument (argc, argv, "-s", skipIdxStr);
    skipIdx = atoi(skipIdxStr.c_str());
  }
  
  shouldRot = find_switch (argc, argv, "-zFix");

  parse_argument (argc, argv, "-cal", hdlCalibration);
  parse_argument (argc, argv, "-pcap", pcapFile);
  
  dirName = pcapFile.substr(0,pcapFile.find_last_of('.',pcapFile.length()));
  mkdir(dirName.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  don = find_switch (argc, argv, "-don");
  cube = find_switch (argc, argv, "-cube");
  intFilt = find_switch (argc, argv, "-intFilt");
  noCloud = find_switch (argc, argv, "-noCloud");
  if (don)
  {
    donDir = dirName + "/don";
    mkdir(donDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  if (cube)
  {
    cubeDir = dirName + "/cube";
    mkdir(cubeDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  if (intFilt)
  {
    parse_argument (argc, argv, "-intFilt", intStr);
    intScale = atof(intStr.c_str());
    intDir = dirName + "/int";
    mkdir(intDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }
  

  MyHdlGrabber grabber (hdlCalibration,pcapFile);
  grabber.start();
  
  
  //skip first cloud, radially starts at a random location
  if (grabber.hasClouds())
  {
    cloud_ = grabber.getNextCloud();
  }

  while (grabber.hasClouds())
  {
    while (cloudIdx < skipIdx)
    {
      cloud_ = grabber.getNextCloud();
      cloudIdx++;
    }
    cloud_ = grabber.getNextCloud();

    if (shouldRot)
    {
      aligned.clear();
      pcl::transformPointCloud (*cloud_,aligned,transMat);
      cloud_->clear();
      (*cloud_) += aligned;
      aligned.clear();
    }

    sprintf(fName,"./%s/%04d.pcd",dirName.c_str(),cloudIdx);
    //save the timestamp to ####.ts
    sprintf(fName2,"./%s/%04d.ts",dirName.c_str(),cloudIdx);
    if (!noCloud)
    {
      cout << "Saving " << fName << " to file...\n";
      pcl::io::savePCDFileBinaryCompressed(fName,*cloud_);
      tStampFile = fopen(fName2,"w");
      if (tStampFile)
      {
	fprintf(tStampFile,"%ld",cloud_->header.stamp);
	fclose(tStampFile);
      }
    }
    
    if (don)
    {
      donFilter.filter(cloud_,donCloud);
      
      sprintf(fName,"./%s/%04d.pcd",donDir.c_str(),cloudIdx);
      cout << "Saving " << fName << " to file...\n";
      pcl::io::savePCDFileBinaryCompressed(fName,*donCloud);
    }
    
    if (cube)
    {
      CubeAligner cubeAlign(cloud_,0.5);
      cubeAlign.identifyClusters();
      std::vector<CubeAligner::CCluster>::iterator clustIt;
      intcloud_->clear();
      int iVal = 0;

      for(clustIt = cubeAlign.getClusters()->begin(); 
        clustIt < cubeAlign.getClusters()->end();
        clustIt++)
      {
        (*intcloud_) += *(*clustIt).toPolygon(iVal++);
      }
      
      if (intcloud_->empty())
      {
	PointXYZI pt;
	pt.x=0;
	pt.y=0;
	pt.z=0;
	intcloud_->push_back(pt);
      }

      sprintf(fName,"./%s/%04d.pcd",cubeDir.c_str(),cloudIdx);
      cout << "Saving " << fName << " to file...\n";
      pcl::io::savePCDFileBinaryCompressed(fName,*intcloud_);
    }
    
    if (intFilt)
    {
      intcloud_->clear();
      donFilter.intensityFilter(intScale,cloud_,intcloud_);
      sprintf(fName,"./%s/%04d.pcd",intDir.c_str(),cloudIdx);
      cout << "Saving " << fName << " to file...\n";
      pcl::io::savePCDFileBinaryCompressed(fName,*intcloud_);
    }
    cloudIdx++;
    
  }
  return (0);
}
