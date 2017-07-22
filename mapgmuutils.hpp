#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/foreach.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <sstream>
#include <iostream>
#include <fstream>
#include <GeographicLib/LambertConformalConic.hpp> 
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace std;
using namespace pcl;
using namespace GeographicLib; 
using namespace Eigen;

#ifndef MAP_GMU_UTILS
#define MAP_GMU_UTILS

namespace mapgmu{

  template<typename PointType>
  class CloudBoundaryFilter
  {
    public:
      CloudBoundaryFilter (
	double x, double X,
	double y, double Y,
	double z, double Z)
      {
	minX = x;
	maxX = X;
	minY = y;
	maxY = Y;
	minZ = z;
	maxZ = Z;
      }
      
      inline int filterCloud(typename PointCloud<PointType>::Ptr cloudIn,
			typename PointCloud<PointType>::Ptr cloudOut)
      {
	typename PointCloud<PointType>::iterator cloudItr;
	
	cloudOut->clear();
	for (cloudItr = cloudIn->begin(); cloudItr != cloudIn->end(); cloudItr++)
	{
	  if ((*cloudItr).x < minX || (*cloudItr).x  > maxX ||
	    (*cloudItr).y < minY || (*cloudItr).y  > maxY ||
	    (*cloudItr).z < minZ || (*cloudItr).z  > maxZ )
	  {
	    continue;
	  }
	  cloudOut->push_back(*cloudItr);
	}
	
	return cloudOut->size();
      }

    private:
      double minX;
      double maxX;
      double minY;
      double maxY;
      double minZ;
      double maxZ;
  };
  
  class PathFinder
  {
    public:
      PathFinder(){}
      
      PathFinder(int nPoints):
	path_vec(),
	segment_lens(),
	currIdx(0),
	ignore_z(false)
      {
	path_vec.resize(nPoints);
	segment_lens.resize(nPoints-1);
      }
      
      
      //supports kml file with LineString (needs to be updated for seg_times to work)
      // or gx:MultiTrack.gx:Track with when times
      PathFinder(string kmlFileName, string lineStringName):
	path_vec(),
	segment_lens(),
	currIdx(0),
	ignore_z(false)
      {
	using boost::property_tree::ptree;
	
	stringstream ss;
	int dint;
	double ddouble;

	// Set up the input datetime format.
	boost::local_time::local_time_input_facet *input_facet 
	    = new boost::local_time::local_time_input_facet("%Y-%m-%dT%H:%M:%sZ");
	ss.imbue(std::locale(ss.getloc(), input_facet));
	boost::local_time::local_date_time ldt(boost::date_time::not_a_date_time);
	boost::posix_time::ptime time_t_epoch(boost::gregorian::date(1970,1,1));
	double epochTime;
	
	vector<string> coord_strs;
	vector<string> point_strs;
	string coordString;
	ptree pt;
	ifstream fileStream;
	
	// Set up basic projection for converting gps to local coordinates
	const double
	a = Constants::WGS84_a<double>(),
	f = 1/298.257222101,                      // GRS80
	lat1 = 40 + 58/60.0, lat2 = 39 + 56/60.0, // standard parallels
	k1 = 1,                                   // scale
	lon0 =-77 - 45/60.0; // origin
	double x,y,z,x0=-1,y0=-1,z0=-1,lat,lon;
	const LambertConformalConic MyConic(a, f, lat1, lat2, k1);
	
	fileStream.open(kmlFileName.c_str(),ios::in);
	if (!fileStream.is_open())
	  throw 10;
	read_xml(fileStream,pt);
	fileStream.close();
	
	boost::optional<ptree&> lsTree;
	ptree trackTree;
	std::pair<ptree::assoc_iterator,ptree::assoc_iterator> pMarks;
	std::pair<ptree::assoc_iterator,ptree::assoc_iterator> gxPts;
	pMarks = pt.get_child("kml.Document").equal_range("Placemark");
	
	for (ptree::const_assoc_iterator pItr(pMarks.first); pItr != pMarks.second; ++pItr)
	{
	  if (pItr->second.get_child_optional("LineString"))
	  {
	    //cout << lsTree->get<string>("name") << "\n";
	    if (pItr->second.get<string>("name") == lineStringName)
	    {
	      coordString = pItr->second.get<string>("LineString.coordinates","").data();
	      boost::split(coord_strs,coordString,boost::is_any_of("\n\t "));
	      
	      for (vector<string>::iterator it=coord_strs.begin();it!=coord_strs.end();it++)
	      {
		boost::split(point_strs,*it,boost::is_any_of(","));
		if (point_strs.size() == 3)
		{
		  lon = atof(point_strs.at(0).c_str());
		  lat = atof(point_strs.at(1).c_str());
		  z = atof(point_strs.at(2).c_str());
		  MyConic.Forward(lon0,lat,lon,x,y);
		  if (x0==-1&&y0==-1&&z0==-1){x0=x;y0=y;z0=z;}
		  x-=x0;y-=y0;z-=z0;
		  //cout << "adding waypoint, x=" << x <<", y=" << y << ", z=" << z << endl;
		  addWayPointXYZ(x,y,z);
		}
		else
		{
		  cout << "wrong number of coordinates: " << point_strs.size() << endl;
		  throw 10;
		}
	      }
	      break;
	    }
	  }
	  else if (pItr->second.get_child_optional("gx:MultiTrack.gx:Track"))
	  {
	    trackTree = pItr->second.get_child("gx:MultiTrack.gx:Track");
	    gxPts = trackTree.equal_range("gx:coord");
	    for (ptree::const_assoc_iterator itr(gxPts.first); itr != gxPts.second; ++itr)
	    {
	      boost::split(point_strs,itr->second.data(),boost::is_any_of(" "));
	      if (point_strs.size() == 3)
	      {
		lon = atof(point_strs.at(0).c_str());
		lat = atof(point_strs.at(1).c_str());
		z = atof(point_strs.at(2).c_str());
		MyConic.Forward(lon0,lat,lon,x,y);
		if (x0==-1&&y0==-1&&z0==-1){x0=x;y0=y;z0=z;}
		x-=x0;y-=y0;z-=z0;
		//cout << "adding waypoint, x=" << x <<", y=" << y << ", z=" << z << endl;
		addWayPointXYZ(x,y,z);
	      }
	      else
	      {
		cout << "wrong number of coordinates: " << point_strs.size() << endl;
		throw 10;
	      }
	    }
	    gxPts = trackTree.equal_range("when");
	    for (ptree::const_assoc_iterator itr(gxPts.first); itr != gxPts.second; ++itr)
	    {
	      ss.str("");
	      ss.clear();
	      ss << itr->second.data();
	      ss >> ldt;
	      
	      boost::posix_time::time_duration diff = ldt.utc_time() - time_t_epoch;
	      epochTime = diff.total_milliseconds()/1000.0;
	      if (seg_times.empty())
	      {
		seg_times.push_back(epochTime);
	      }
	      else
	      {
		seg_times.push_back(epochTime - seg_times.at(0));
	      }
	    }
	    if (!seg_times.empty())
	      seg_times.at(0) = 0;
	  }
	}
      }
      
      //load the kml after constructor
      inline void loadKml(string kmlFileName, string lineStringName)
      {
	unsigned int i;
	//load with constructor
	PathFinder pf(kmlFileName,lineStringName);
	
	//deep copy everything out of pf
	currIdx = 0;
	segment_lens.resize(0);
	path_vec.resize(0);
	for (i=0;i<pf.path_vec.size();i++)
	{
	  addWayPoint(pf.path_vec.at(i));
	}
	for (i=0;i<pf.seg_times.size();i++)
	{
	  seg_times.push_back(pf.seg_times.at(i));
	}
      }
      
      inline void setIgnoreZ(){ignore_z=true;}
      
      inline void addWayPoint(PointXYZ point)
      {
	Eigen::Vector3f diff;
	path_vec.push_back(point);
	if (currIdx > 0)
	{
	  diff = path_vec.at(currIdx).getVector3fMap()-
		path_vec.at(currIdx-1).getVector3fMap();
	  if (ignore_z){diff(2)=0;}
	  segment_lens.push_back(sqrt(diff.norm()));
	}
	currIdx++;
      }
      
      inline void addWayPointXY(double x, double y)
      {
	addWayPointXYZ(x,y,0.0);
      }
      
      inline void addWayPointXYZ(double x, double y, double z)
      {
	PointXYZ newPoint;
	newPoint.x=x;
	newPoint.y=y;
	newPoint.z=z;
	
	addWayPoint(newPoint);
      }
      
      inline double distanceFromWaypoint(int idx, PointXYZ p1)
      {
	Eigen::Vector3f diff = path_vec.at(idx).getVector3fMap()-
		p1.getVector3fMap();
	if (ignore_z)
	{
	  diff(2) = 0;
	}
	return sqrt(diff.norm());
      }
      
      inline bool isInSegment(int idx, PointXYZ p1)
      {
	if (idx > currIdx)
	  return false;
	else
	  return distanceFromWaypoint(idx,p1)<=segment_lens.at(idx);
      }
      
      //produces the vector transform that rotates the vetor from p1->p2
      // onto the segment of the track denoted by segIdx
      inline Matrix4f rotateAlignment(PointXYZ p1, PointXYZ p2, int segIdx)
      {
	Matrix4f retMat = Matrix4f::Identity();
	Matrix3f rotMat;
	int i,j;
	if (segIdx > (int)path_vec.size()  + 1 )
	{
	  cout << "Path has no completed segment #" << segIdx << endl;
	  throw 10;
	}
	Vector3f myVec = (p2.getVector3fMap() - p1.getVector3fMap());
	if (ignore_z){myVec(2) = 0;}
	Vector3f segVec = (path_vec.at(segIdx+1).getVector3fMap() - 
	    path_vec.at(segIdx).getVector3fMap());
	Vector3f crossVec = myVec.cross(segVec);
	crossVec /= crossVec.norm();
	double aCosVal = myVec.dot(segVec)/(myVec.norm()*segVec.norm());
	double theta = acos(aCosVal);
	
	cout << "rotating " << p1 << "->"<<p2<<endl;
	cout << "onto " << path_vec.at(segIdx) << "->"<<path_vec.at(segIdx+1)<<endl;
	cout << "acos = " << aCosVal << ", theta= " << theta << endl;
	
	
	rotMat = AngleAxisf(theta,crossVec);
	for(i=0;i<3;i++){for(j=0;j<3;j++){
	  retMat(i,j) = rotMat(i,j);
	}}
	
	Vector3f norm1,norm2;
	norm1 = segVec/segVec.norm();
	norm2 = rotMat*myVec;
	norm2 /= norm2.norm();
	
	//cout << "Normalized segment = " << norm1 << endl;
	//cout << "Normalized rotated path = " << norm2 << endl;
	
	//cout << "align mat: " << retMat << endl;
	
	return retMat;
      }
      
      //same as rotateAlignment, doesn't use Z
      inline Matrix4f rotateAlignment2D(PointXYZ p1, PointXYZ p2, int segIdx)
      {
	Matrix4f retMat = Matrix4f::Identity();
	Matrix3f rotMat;
	int i,j;
	if (segIdx > (int)path_vec.size()  + 1 )
	{
	  cout << "Path has no completed segment #" << segIdx << endl;
	  throw 10;
	}
	Vector3f myVec = (p2.getVector3fMap() - p1.getVector3fMap());
	myVec(2) = 0;
	Vector3f segVec = (path_vec.at(segIdx+1).getVector3fMap() - 
	    path_vec.at(segIdx).getVector3fMap());
	segVec(2) = 0;
	Vector3f crossVec = myVec.cross(segVec);
	crossVec /= crossVec.norm();
	double aCosVal = myVec.dot(segVec)/(myVec.norm()*segVec.norm());
	double theta = acos(aCosVal);
	
	rotMat = AngleAxisf(theta,crossVec);
	for(i=0;i<3;i++){for(j=0;j<3;j++){
	  retMat(i,j) = rotMat(i,j);
	}}
	
	Vector3f norm1,norm2;
	norm1 = segVec/segVec.norm();
	norm2 = rotMat*myVec;
	norm2 /= norm2.norm();
	
	//cout << "Normalized segment = " << norm1 << endl;
	//cout << "Normalized rotated path = " << norm2 << endl;
	
	//cout << "align mat: " << retMat << endl;
	
	return retMat;
      }
      
      inline PointXYZ extractOriginD(Eigen::Matrix4d rotMat)
      {
	PointXYZ newPoint;
	newPoint.x = rotMat(0,3);
	newPoint.y = rotMat(1,3);
	newPoint.z = rotMat(2,3);
	
	
	return newPoint;
      }
      
      inline PointXYZ extractOrigin(Eigen::Matrix4f rotMat)
      {
	PointXYZ newPoint;
	newPoint.x = rotMat(0,3);
	newPoint.y = rotMat(1,3);
	newPoint.z = rotMat(2,3);
	
	
	return newPoint;
      }
      
      //creates a scaled vector for time extending dur after
      // time sec, does not account for wrapping past end of segments
      // time
      inline Vector3f getDeltaForTime(double seconds,double dur)
      {
	int i;
	double delta;
	PointXYZ p1,p2;
	Vector3f retMat;
	
    cout<<"seg_times.size="<<seg_times.size()<<endl;

	for (i=1;i<seg_times.size();i++)
	{
	  if (seconds < seg_times.at(i))
	  {
	    break;
	  }
	}
	if (i>=seg_times.size())
	{
	  cout << "Path doesn't have points at time " << seconds << endl;
	  throw 10;
	}
	else
	{
	  delta = (seg_times.at(i)-seg_times.at(i-1));
	  p1 = path_vec.at(i-1);
	  p2 = path_vec.at(i);
	  retMat(0) = ((p2.x-p1.x)/delta)*dur;
	  retMat(1) = ((p2.y-p1.y)/delta)*dur;
	  retMat(2) = ((p2.z-p1.z)/delta)*dur;
	  
	  return retMat;
	}
      }
      
      //finds the segment that seconds falls within and the point
      // within the segment (assuming constant vel) 
      // and creates a transformation matrix that rotates (0,1,0)
      // onto the segment and translates the distance of the vector
      // from the start of the segment to the point found inside the segment
      //Used in SimplePCDViewer to cause the single frames to follow the
      // track of the GPS, assuming that the gps is accurate, (0,1,0) is the
      // forward motion of the cart
      inline Matrix4f getTransForTime(double seconds)
      {
	int i;
	double delta,delta2;
	Matrix4f retMat;
	PointXYZ p1,p2;
	
	p1.x=0;p1.y=0;p1.z=0;
	p2.x=0;p2.y=1;p2.z=0;
	
	for (i=1;i<seg_times.size();i++)
	{
	  if (seconds < seg_times.at(i))
	  {
	    break;
	  }
	}
	if (i>=seg_times.size())
	{
	  cout << "Path doesn't have points at time " << seconds << endl;
	  throw 10;
	}
	else
	{
	  retMat = rotateAlignment2D(p1,p2,i-1);
	  delta = (seg_times.at(i)-seg_times.at(i-1));
	  delta2 = (seconds-seg_times.at(i-1));
	  p1 = path_vec.at(i-1);
	  p2 = path_vec.at(i);
	  retMat(0,3) = p1.x + ((p2.x-p1.x)/delta)*delta2;
	  retMat(1,3) = p1.y + ((p2.y-p1.y)/delta)*delta2;
	  retMat(2,3) = p1.z + ((p2.z-p1.z)/delta)*delta2;
	  
	  return retMat;
	}
      }
      
      inline PointXYZ getWaypoint(int idx)
      {
	if (idx >= (int)path_vec.size() )
	{
	  cout << "Path has no completed segment #" << idx << endl;
	  throw 10;
	}
	return path_vec.at(idx);
      }
      
      inline int getWaypointCount(void)
      {
	return (int) path_vec.size();
      }
      
    private:
      vector<PointXYZ> path_vec;
      vector<double> segment_lens;
      vector<double> seg_times;
      double deltaRad;
      int currIdx;
      bool ignore_z;
  };
}

#endif
