#include <pcl/search/kdtree.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>
#include <Eigen/Dense>
#include <cmath>

#ifndef DON_ESTIMATOR_HPP
#define DON_ESTIMATOR_HPP

#ifndef MY_PI
#define MY_PI 3.14159265
#endif

using namespace std;
using namespace pcl;
using namespace Eigen;


//class used to calculate the normals on a cloud of points
// DoN Filtering done by DoNFilter

template<typename PointType>
class DonEstimator
{
public:
  typedef PointCloud<PointType> Cloud;
  typedef typename PointCloud<PointType>::Ptr CloudPtr;
  typedef typename PointCloud<PointType>::iterator CloudItr;
  typedef PointNormal PNorm;
  typedef PointCloud<PNorm> NormalCloud;
  typedef typename PointCloud<PNorm>::iterator NCloudItr;
  typedef typename PointCloud<PNorm>::Ptr NormalCloudPtr;
  
  DonEstimator()
  : maxNeighbors_(100)
  , minDistThresh_(sin((28.6/64)*MY_PI/180.0)*3.5)
  , maxDistScale_(3.5)
  , dotProdThresh_(cos(45.0*MY_PI/180.0))
  , viewpoint_(0.0,0.0,0.0)
  , myMesh_(0)
  , createMesh(false){}
  
  inline void setInputCloud(CloudPtr cloud){inputCloud_ = cloud;}
  
  inline void setMinDistThread(double dist){minDistThresh_=dist;}
  
  inline void setMaxDistScale(double scale){maxDistScale_=scale;}
  
  inline void setMaxNeigbors(int max){maxNeighbors_=max;}
  
  inline void setDotProductThreshold(double value){dotProdThresh_=value;}
  
  inline void setViewpoint(float x, float y, float z)
  {
    viewpoint_(0) = x;
    viewpoint_(1) = y;
    viewpoint_(2) = z;
  }
  
  inline std::vector<pcl::Vertices> buildMesh(CloudPtr baseCloud)

  {
    NormalCloudPtr dummy(new NormalCloud);
    setInputCloud(baseCloud);
    calculateDoN(dummy);
    
    return myMesh_;
  }
  
  inline double getMaxRadiusForPoint(PointType pt)
  {
    Vector3f distToVp;
    double minDist;
    
    distToVp = pt.getVector3fMap()-viewpoint_;
    minDist = sqrt(distToVp.norm())*minDistThresh_;
    
    return minDist*maxDistScale_;
  }
  
  inline double getMinRadiusForPoint(PointType pt)
  {
    Vector3f distToVp;
    double minDist;
    
    distToVp = pt.getVector3fMap()-viewpoint_;
    minDist = sqrt(distToVp.norm())*minDistThresh_;
    
    return minDist;
  }
  
  inline void calculateDoN(NormalCloudPtr cloudOut)
  {
    CloudItr cItr;
    NCloudItr nItr;
    int i,res,rmn = (maxNeighbors_+1);
    PNorm smallNorm,dNorm;
    vector<int> ptIds(rmn);
    vector<float> nDists(rmn);
    Vector3f distToVp;
    Vector4f nBuff;
    double minDist,maxDist;
    int badCount = 0;
    
    cloudOut->clear();
    
    myTree_.setInputCloud(inputCloud_);
    
    for (cItr=inputCloud_->begin();cItr<inputCloud_->end();cItr++)
    {
      //reset diff and p2->p5
      dNorm.curvature = 0;
      dNorm.x = cItr->x;
      dNorm.y = cItr->y;
      dNorm.z = cItr->z;
      
      
      //calculate distance to viewpoint
      distToVp = cItr->getVector3fMap()-viewpoint_;
      minDist = sqrt(distToVp.norm())*minDistThresh_;
      maxDist = minDist*maxDistScale_;
      minDist *= minDist;
      
      //do nearest neighbor search
      res = myTree_.radiusSearch(*cItr,maxDist,ptIds,nDists,rmn);
      if (!res){cout << "Found 0 Points!!!!\n";}
      for ( i=1; i<res; i++ )
      {
	if (nDists[i] > minDist)
	{
	  break;
	}
      }
      
      vector<int> normSet(ptIds.begin(),ptIds.begin() + res);
      
      computePointNormal(*inputCloud_,normSet,nBuff,dNorm.curvature);
      flipNormalTowardsViewpoint(*cItr,viewpoint_[0],viewpoint_[1],viewpoint_[2],nBuff);
      dNorm.data_n[0] = nBuff[0];
      dNorm.data_n[1] = nBuff[1];
      dNorm.data_n[2] = nBuff[2];
      
      normSet.resize(i);
      
      computePointNormal(*inputCloud_,normSet,nBuff,smallNorm.curvature);
      flipNormalTowardsViewpoint(*cItr,viewpoint_[0],viewpoint_[1],viewpoint_[2],nBuff);
      smallNorm.data_n[0] = nBuff[0];
      smallNorm.data_n[1] = nBuff[1];
      smallNorm.data_n[2] = nBuff[2];
      
      dNorm.getNormalVector3fMap() = (smallNorm.getNormalVector3fMap() - dNorm.getNormalVector3fMap())/2.0;
      
      if( !pcl_isfinite (dNorm.normal_x) ||
	  !pcl_isfinite (dNorm.normal_y) ||
	  !pcl_isfinite (dNorm.normal_z) )
      {
	badCount++;
      }
      else
      {
	dNorm.curvature = dNorm.getNormalVector3fMap ().norm();
	cloudOut->push_back(dNorm);
      }
    }
    
    //cout << "got " << badCount << " bad points out of " << inputCloud_->size() << endl;
  }
  
  inline void setMeshFlag(bool flag){createMesh = flag;}
    
protected:
  
  KdTreeFLANN<PointType> myTree_;
  KdTreeFLANN<PNorm> nTree_;
  CloudPtr inputCloud_;
  int maxNeighbors_;
  double minDistThresh_;
  double maxDistScale_;
  double dotProdThresh_;
  Vector3f viewpoint_;
  vector<Vertices> myMesh_;
  bool createMesh;
};

#endif