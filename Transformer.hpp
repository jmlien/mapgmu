#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <string>
#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <Eigen/Dense>
#include <pcl/registration/icp.h>
#include <pcl/registration/ndt.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/registration/icp_nl.h>
#include <pcl/search/organized.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/features/fpfh.h>
#include "DoNFilter.hpp"
#include "CubeAligner.hpp"
#include "mapgmuutils.hpp"

#ifndef TRANSFORMER_H
#define TRANSFORMER_H

using namespace pcl;
using namespace std;

typedef pcl::PointXYZI PointType;
typedef PointCloud<pcl::PointXYZI> Cloud;
typedef PointCloud<pcl::PointXYZI>::Ptr CloudPtr;
typedef pcl::PointNormal PointNormalT;
typedef pcl::PointCloud<PointNormalT> PointCloudWithNormals;
typedef pcl::PointCloud<PointNormalT>::Ptr PointCloudWithNormalsPtr;
typedef pcl::KdTreeFLANN<PointType> TreeType;
typedef TreeType::Ptr TreePtr;

class MyPointRepresentation : public pcl::PointRepresentation <PointNormalT>
{
  using pcl::PointRepresentation<PointNormalT>::nr_dimensions_;
public:
  MyPointRepresentation ()
  {
    // Define the number of dimensions
    nr_dimensions_ = 4;
  }

  // Override the copyToFloatArray method to define our feature vector
  virtual inline void copyToFloatArray (const PointNormalT &p, float * out) const
  {
    // < x, y, z, curvature >
    out[0] = p.x;
    out[1] = p.y;
    out[2] = p.z;
    out[3] = p.curvature;
  }
};

class Transformer
{
public:
  
  inline double getFitness(CloudPtr base, CloudPtr toAlign, Eigen::Matrix4f mat)
  {
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    
    icp.setMaxCorrespondenceDistance (0.5);
    icp.setMaximumIterations (0);
    
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final);
    return icp.getFitnessScore();
  }
  
  inline Eigen::Matrix4f getCubeTransform(CloudPtr base, CloudPtr toAlign,
					  double cSize, double *eVal)
  {
    Eigen::Matrix4f baseTrans = Eigen::Matrix4f::Identity();
    CloudPtr copyCloud(new Cloud);
    CubeAligner cb(base,cSize);
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    return cb.alignToMatrix(baseTrans,copyCloud,eVal);
  }
  
  inline Eigen::Matrix4f getCubeICPTransform(CloudPtr base, CloudPtr toAlign, 
					     double cSize, double *eVal)
  {
    Eigen::Matrix4f tmpTrans;
    Eigen::Matrix4f baseTrans = Eigen::Matrix4f::Identity();
    
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    CloudPtr copyCloud(new Cloud);
    int maxIter = 100;
    double eps = 1e-8;
    
    CubeAligner cb(base,cSize);
    
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    //binary cube alignment using cubic filtered toAlign
    tmpTrans =  cb.alignToMatrix(baseTrans,copyCloud,eVal);
    
    //scale added for tuning, not necessary
    icp.setMaxCorrespondenceDistance (cSize*2.0);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    
    //icp using the filtered copyCloud
    icp.setInputSource(copyCloud);
    icp.setInputTarget(base);
    icp.align(Final,tmpTrans);
    
    return icp.getFinalTransformation();
  }
  
  //assumes gpsGuess already identity, sets the translation values
  inline void calcGpsGuess(Eigen::Matrix4f& lastTrans,mapgmu::PathFinder *pf,
			   double currTime, double dt, Eigen::Matrix4f& gpsGuess)
  {
    Vector3f rotOffset;
    //create gps estimate, store in 4x4 matrix
    rotOffset =  lastTrans.block(0,0,3,3).transpose()*pf->getDeltaForTime(currTime,dt);
    gpsGuess(0,3)=rotOffset(0);
    gpsGuess(1,3)=rotOffset(1);
    gpsGuess(2,3)=rotOffset(2);
  }
  
  inline Eigen::Matrix4f getGpsIcpNoCubeTransform(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, double *cSize,
					  double currTime, double dt)
  {    
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp2;
    pcl::PointCloud<pcl::PointXYZI> Final;
    int maxIter = 50;
    double eps = 1e-8;
    double eVal = *cSize*2.0;
    
    Matrix4f gpsGuess = Matrix4f::Identity();
    Matrix4f transMat;
    
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    //icp pass 1
    icp.setMaxCorrespondenceDistance (eVal);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final,gpsGuess);
    
    transMat =  icp.getFinalTransformation();
    
    Final.clear();
    
    //icp pass 2, smaller corr, higher euclid fitness
    icp2.setMaxCorrespondenceDistance (*cSize/15.0);
    icp2.setMaximumIterations (maxIter);
    icp2.setEuclideanFitnessEpsilon (eps*10);
    icp2.setInputSource(toAlign);
    icp2.setInputTarget(base);
    icp2.align(Final,transMat);
    
    transMat = icp2.getFinalTransformation();
    
    return transMat;
  }
  
  inline Eigen::Matrix4f getGpsIcpTransform(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, double *cSize,
					  double currTime, double dt)
  {    
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp2;
    pcl::PointCloud<pcl::PointXYZI> Final;
    CloudPtr copyCloud(new Cloud);
    int maxIter = 50;
    double eps = 1e-8;
    double eVal = *cSize*2.0;
    
    Vector3f rotOffset;
    Matrix4f gpsGuess = Matrix4f::Identity();
    Matrix4f transMat;
    CubeAligner cb(base,*cSize);
    
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    CVertex c1(gpsGuess);
    
    //first pass icp, filtered copyCloud used
    icp.setMaxCorrespondenceDistance (eVal);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(copyCloud);
    icp.setInputTarget(base);
    icp.align(Final,gpsGuess);
    
    transMat =  icp.getFinalTransformation();
    
    Final.clear();
    
    //second pass icp, full toAlign cloud used, smalled corr, higher euc
    icp2.setMaxCorrespondenceDistance (*cSize/15.0);
    icp2.setMaximumIterations (maxIter);
    icp2.setEuclideanFitnessEpsilon (eps*10);
    icp2.setInputSource(toAlign);
    icp2.setInputTarget(base);
    icp2.align(Final,transMat);
    
    transMat = icp2.getFinalTransformation();
    
    //uncomment for distance/metric calculations
    /*CVertex c2(transMat);
    CVertex c3 = c2-c1;
    
    std::cout << "Cube: " << *cSize << std::endl 
		<< "Conv: " << eVal << std::endl
		<< "Converged: " << icp2.hasConverged () << std::endl
		<< "Fitness: " << icp2.getFitnessScore () << std::endl
		<< "Vexp: " << c1.getDist() << std::endl
		<< "Vfit: " << c2.getDist() << std::endl
		<< "PathDist: " << c3.getDist()  << std::endl;*/
    
    return transMat;
  }
  
  inline Eigen::Matrix4f getGpsIcpCubeTransform(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, double *cSize,
					  double *eVal, double currTime, double dt)
  {    
    Eigen::Matrix4f tmpTrans;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    CloudPtr copyCloud(new Cloud);
    int maxIter = 50;
    double eps = 1e-6;
    
    Vector3f rotOffset;
    Matrix4f gpsGuess = Matrix4f::Identity();
    CubeAligner cb(base,*cSize);
    
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    //icp pass
    icp.setMaxCorrespondenceDistance (*cSize);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(copyCloud);
    icp.setInputTarget(base);
    icp.align(Final,gpsGuess);
    
    std::cout << "Converged: " << icp.hasConverged () << std::endl
		<< "Fitness: " << icp.getFitnessScore () << std::endl;
    
    tmpTrans = icp.getFinalTransformation();
    
    //cubic alignment
    return cb.alignToMatrix(tmpTrans,toAlign,eVal);
  }
  
  inline Eigen::Matrix4f getCubeGpsIcpTransform(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, double *cSize,
					  double *eVal, double currTime, double dt)
  {
    Eigen::Matrix4f tmpTrans;
    
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp2;
    pcl::PointCloud<pcl::PointXYZI> Final;
    CloudPtr copyCloud(new Cloud);
    int maxIter = 100;
    double eps = 1e-8;
    
    Vector3f rotOffset;
    Matrix4f gpsGuess = Matrix4f::Identity();
    Matrix4f transMat;
    CubeAligner cb(base,*cSize);
    
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    tmpTrans =  cb.alignToMatrix(gpsGuess,copyCloud,eVal);
    
    //icp pass one, filtered copyCloud
    icp.setMaxCorrespondenceDistance (*cSize*2.0);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(copyCloud);
    icp.setInputTarget(base);
    icp.align(Final,tmpTrans);
    
    transMat = icp.getFinalTransformation();
    
    Final.clear();
    
    //icp pass two full toAlign cloud
    icp2.setMaxCorrespondenceDistance (*cSize/15.0);
    icp2.setMaximumIterations (maxIter);
    icp2.setEuclideanFitnessEpsilon (eps*100);
    icp2.setInputSource(toAlign);
    icp2.setInputTarget(base);
    icp2.align(Final,transMat);
    
    transMat = icp2.getFinalTransformation();
    
    //attempts at dynamically updating eVal
    /* *eVal = (*eVal + cb.getNextEval(&gpsGuess,&tmpTrans))/2.0;
    if (*eVal < 0.05) *eVal = 0.05;
    
    std::cout << "New EVAL " << *eVal << std::endl;
    
    std::cout << "Converged: " << icp.hasConverged () << std::endl
		<< "Fitness: " << icp.getFitnessScore () << std::endl;*/
    
    return transMat;
  }
  
  inline Eigen::Matrix4f getCubeGpsIcpTransformAdapt(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, CVertex lastVel,
					  double currTime, double dt)
  {
    Eigen::Matrix4f cubeTrans;
    Eigen::Vector3f vel = pf->getDeltaForTime(currTime,dt);
    
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp2;
    pcl::PointCloud<pcl::PointXYZI> Final;
    CloudPtr copyCloud(new Cloud);
    int maxIter = 100;
    double eps = 1e-8;
    
    Vector3f rotOffset;
    
    Matrix4f transMat;
    double cSize;
    double eVal;
    
    Matrix4f gpsGuess = Matrix4f::Identity();
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    cSize = vel.norm();
    if (lastVel.getDist()>0)
    {
      cSize += lastVel.getDist();
      cSize /= 2.0;
    }
    
    if (cSize > 0.45)cSize = 0.45;
    if (cSize < 0.05)cSize = 0.05;
    
    eVal = cSize / 2.0;
    
    CubeAligner cb(base,cSize);
    
    *copyCloud += *toAlign;
    cb.filterCloudCubic(copyCloud);
    
    //cube alignment with gps guess seed
    cubeTrans =  cb.alignToMatrix(gpsGuess,copyCloud,&eVal);
    
    //icp pass one with cube alignment seed, filtered copyCloud
    icp.setMaxCorrespondenceDistance (cSize*1.33);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(copyCloud);
    icp.setInputTarget(base);
    icp.align(Final,cubeTrans);
    
    transMat = icp.getFinalTransformation();
    
    Final.clear();
    
    //icp pass two full toAlign cloud
    icp2.setMaxCorrespondenceDistance (cSize/10.0);
    icp2.setMaximumIterations (maxIter/5.0);
    icp2.setEuclideanFitnessEpsilon (eps*100);
    icp2.setInputSource(toAlign);
    icp2.setInputTarget(base);
    icp2.align(Final,transMat);
    
    transMat = icp2.getFinalTransformation();
    
    return transMat;
  }
    
  inline Eigen::Matrix4f getCubeGpsTransform(CloudPtr base, CloudPtr toAlign,
					  Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf, double *cSize,
					  double *eVal, double currTime, double dt)
  {        
    CubeAligner cb(base,*cSize);
    
    Matrix4f gpsGuess = Matrix4f::Identity();
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    CVertex rv(gpsGuess);
    *eVal = rv.getDist();
    
    return cb.alignToMatrix(gpsGuess,toAlign,eVal);
  }
  
    
  inline Eigen::Matrix4f getGpsTransform(Eigen::Matrix4f trans,
					  mapgmu::PathFinder *pf,
					  double currTime, double dt)
  {        
    Matrix4f gpsGuess = Matrix4f::Identity();
    calcGpsGuess(trans,pf,currTime,dt,gpsGuess);
    
    return gpsGuess;
  }
  
  inline Eigen::Matrix4f getICPTransform(CloudPtr base, CloudPtr toAlign)
  {
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    int maxIter = 50;
    double eps = 1e-6;
    double eVal = 1.0;
    
    icp.setMaxCorrespondenceDistance (eVal);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final);
    
    return icp.getFinalTransformation();
  }
  
  inline Eigen::Matrix4f getICPVoxTransform(CloudPtr base, CloudPtr toAlign)
  {
    pcl::IterativeClosestPoint<PointType, PointType> icp;
    pcl::VoxelGrid<PointType> voxGrid;
    CloudPtr voxBase (new Cloud);
    CloudPtr voxAlign (new Cloud);
    Cloud Final;
    double maxCorrDist = 0.2;
    int maxIter = 50;
    double eps = 1e-8;
    double leafSize = 0.1;
    
    char *corr,*epsC, *leafStr;
    
    corr = getenv("CORR_DIST");
    epsC = getenv("CONV_EPS");
    leafStr = getenv("LEAF_SIZE");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    if (leafStr)
    {
      leafSize = atof(leafStr);
    }
    
    voxGrid.setLeafSize(leafSize,leafSize,leafSize);
    voxGrid.setInputCloud(base);
    voxGrid.filter(*voxBase);
    voxGrid.setInputCloud(toAlign);
    voxGrid.filter(*voxAlign);
    
    //icp pass, both source and dest voxel filtered
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(voxAlign);
    icp.setInputTarget(voxBase);
    icp.align(Final);
    
    std::cout << maxCorrDist << "," << eps << "," << icp.getFitnessScore () << "," << icp.hasConverged () << ",";
    
    return icp.getFinalTransformation();
  }
  
  // used to determine next single step ICP where [last] is the last alignment
  inline Eigen::Matrix4f getNextICPStep(CloudPtr base, CloudPtr toAlign, 
					Eigen::Matrix4f last)
  {
    pcl::IterativeClosestPoint<PointType, PointType> icp;
    pcl::VoxelGrid<PointType> voxGrid;
    Cloud Final;
    double maxCorrDist = 0.3;
    double eps = 1e-8;
    
    char *corr,*epsC;
    
    corr = getenv("CORR_DIST");
    epsC = getenv("CONV_EPS");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (1);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final,last);
    
    std::cout << "ICP Converged:" << icp.hasConverged () << std::endl
	<< "Fitness: " << icp.getFitnessScore () << std::endl
	<< "EucFit Eps: " << icp.getEuclideanFitnessEpsilon () << std::endl;
    
    return icp.getFinalTransformation();
  }
  
  
  //regular icp transform using "guess" to seed the transform
  inline Eigen::Matrix4f getICPTransformGuess(CloudPtr base, CloudPtr toAlign,
					      Eigen::Matrix4f guess,
					      double maxCorrDist)
  {
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    int maxIter = 50;
    double eps = 1e-6;
    
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (maxIter);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final, guess);
    
    std::cout << "Converged: " << icp.hasConverged () << std::endl
		<< "Fitness: " << icp.getFitnessScore () << std::endl;
    
    return icp.getFinalTransformation();
  }
  
  inline Eigen::Matrix4f getICPIntTransform(CloudPtr base, CloudPtr toAlign)
  {
    pcl::IterativeClosestPoint<PointType, PointType> icp;
    CloudPtr filter_base(new Cloud);
    CloudPtr filter_toAlign(new Cloud);
    Cloud Final;
    Eigen::Matrix4f intensityMat = Eigen::Matrix4f::Identity();
    double maxCorrDist = 2.0;
    double eps = 1e-8;
    DoNFilter<PointType> donFilter;
    
    char *corr,*epsC;
    
    corr = getenv("CORR_DIST");
    epsC = getenv("CONV_EPS");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    
    //filter by intensity, see DoNFilter.hpp
    donFilter.intensityFilter(1.2,base,filter_base);
    donFilter.intensityFilter(1.2,toAlign,filter_toAlign);
    
    //icp pass 1, use filtered clouds
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (50);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(filter_toAlign);
    icp.setInputTarget(filter_base);
    icp.align(Final);
    
    intensityMat = icp.getFinalTransformation();
    
    //icp pass 2, use full clouds
    icp.setMaxCorrespondenceDistance (0.1);
    icp.setInputSource(toAlign);
    icp.setInputTarget(base);
    icp.align(Final,intensityMat);
    
    std::cout << maxCorrDist << "," << eps << "," << icp.getFitnessScore () << "," << icp.hasConverged () << ",";
    
    return icp.getFinalTransformation();
  }
  
  //this was from the PCL site, may work if tweaked, never bother to fix it
  /*inline Eigen::Matrix4f getICPFPFHTransform(CloudPtr base, CloudPtr toAlign,
						      Eigen::Matrix4f guess)
  {
    pcl::IterativeClosestPoint<PointNormalT, PointNormalT> icp;
    PointCloudWithNormalsPtr filter_base(new PointCloudWithNormals);
    PointCloudWithNormalsPtr filter_toAlign(new PointCloudWithNormals);
    pcl::PointCloud<PointNormalT> Final;
    
    pcl::NormalEstimation<PointType, PointNormalT> norm_est;
    pcl::search::KdTree<PointType>::Ptr tree (new pcl::search::KdTree<PointType> ());
    norm_est.setSearchMethod (tree);
    norm_est.setRadiusSearch (0.1);
    
    norm_est.setInputCloud (base);
    norm_est.compute (*filter_base);
    pcl::copyPointCloud (*base, *filter_base);
    
    norm_est.setInputCloud (toAlign);
    norm_est.compute (*filter_toAlign);
    pcl::copyPointCloud (*toAlign, *filter_toAlign);
    
    pcl::FPFHEstimation<PointType, pcl::Normal, pcl::FPFHSignature33> fpfh;
    // alternatively, if cloud is of tpe PointNormal, do fpfh.setInputNormals (cloud);

    fpfh.setSearchMethod (tree);

    // Output datasets
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_base (new pcl::PointCloud<pcl::FPFHSignature33> ());
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhs_toAlign (new pcl::PointCloud<pcl::FPFHSignature33> ());

    // Use all neighbors in a sphere of radius 5cm
    // IMPORTANT: the radius used here has to be larger than the radius used to estimate the surface normals!!!
    fpfh.setRadiusSearch (0.2);

    // Compute the features
    //fpfh.setInputCloud (base);
    fpfh.setInputNormals (filter_base);
    fpfh.compute (*fpfhs_base);
    //fpfh.setInputCloud (toAlign);
    fpfh.setInputNormals (filter_toAlign);
    fpfh.compute (*fpfhs_toAlign);
    
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    // Set the maximum number of iterations (criterion 1)
    icp.setMaximumIterations (maxIter);
    // Set the transformation epsilon (criterion 2)
    icp.setEuclideanFitnessEpsilon (eps);
    
    icp.setInputSource(filter_toAlign);
    icp.setInputTarget(filter_base);
    icp.align(Final,guess);
    
    std::cout << "Converged: " << icp.hasConverged () << std::endl
		<< "Fitness: " << icp.getFitnessScore () << std::endl;
    
    return icp.getFinalTransformation();
    return Eigen::Matrix4f::Identity();
  }*/
  
  inline Eigen::Matrix4f getICPDoNTransform(CloudPtr base, CloudPtr toAlign)
  {
    pcl::IterativeClosestPoint<PointNormalT, PointNormalT> icp;
    PointCloudWithNormalsPtr filter_base(new PointCloudWithNormals);
    PointCloudWithNormalsPtr filter_toAlign(new PointCloudWithNormals);
    pcl::PointCloud<PointNormalT> Final;
    double maxCorrDist = 2.0;
    double eps = 1e-8;
    DoNFilter<PointType> donFilter;
    
    char *corr,*epsC;
    
    corr = getenv("CORR_DIST");
    epsC = getenv("CONV_EPS");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    
    donFilter.filter(base,filter_base);
    donFilter.filter(toAlign,filter_toAlign);
    
    //icp with both source and dest don filtered
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (150);
    icp.setEuclideanFitnessEpsilon (eps);
    icp.setInputSource(filter_toAlign);
    icp.setInputTarget(filter_base);
    icp.align(Final);
    
    std::cout << maxCorrDist << "," << eps << "," << icp.getFitnessScore () << "," << icp.hasConverged () << ",";	
    
    return icp.getFinalTransformation();
  }
  
  //one step DoN, not efficient since DoN filtering is done on each step
  inline Eigen::Matrix4f getNextDoNStep(CloudPtr base, CloudPtr toAlign, 
					Eigen::Matrix4f last)
  {
    pcl::IterativeClosestPoint<PointNormalT, PointNormalT> icp;
    PointCloudWithNormalsPtr filter_base(new PointCloudWithNormals);
    PointCloudWithNormalsPtr filter_toAlign(new PointCloudWithNormals);
    pcl::PointCloud<PointNormalT> Final;
    double maxCorrDist = 0.3;
    int maxIter = 1;
    double eps = 1e-8;
    DoNFilter<PointType> donFilter;
    
    char *corr,*epsC;
    
    corr = getenv("CORR_DIST");
    epsC = getenv("CONV_EPS");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    
    donFilter.filter(base,filter_base);
    donFilter.filter(toAlign,filter_toAlign);
    
    icp.setMaxCorrespondenceDistance (maxCorrDist);
    icp.setMaximumIterations (maxIter);
    icp.setTransformationEpsilon (eps);
    icp.setInputSource(filter_toAlign);
    icp.setInputTarget(filter_base);
    icp.align(Final,last);
    
    std::cout << "ICP DoN Converged:" << icp.hasConverged () << std::endl
	<< "Fitness: " << icp.getFitnessScore () << std::endl
	<< "EucFit Eps: " << icp.getEuclideanFitnessEpsilon () << std::endl;
    
    return icp.getFinalTransformation();
  }
  
  inline Eigen::Matrix4f getNDTDoNTransform(CloudPtr base, CloudPtr toAlign)
  {
    PointCloudWithNormalsPtr filter_base(new PointCloudWithNormals);
    PointCloudWithNormalsPtr filter_toAlign(new PointCloudWithNormals);
    pcl::PointCloud<PointNormalT> Final;
    
    DoNFilter<PointType> donFilter;
    
    char *stepSzC,*iter,*epsC,*resC;
    double stepSize = 0.1;
    int maxIter = 35;
    double eps = 0.01;
    double resolution = 1.0;
    
    stepSzC = getenv("STEP_SIZE");
    iter = getenv("MAX_ITER");
    epsC = getenv("CONV_EPS");
    resC = getenv("NDT_RES");
    if (resC)
    {
      resolution = atof(resC);
    }
    if (stepSzC)
    {
      stepSize = atof(stepSzC);
    }
    if (iter)
    {
      maxIter =atoi(iter);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    
    donFilter.filter(base,filter_base);
    donFilter.filter(toAlign,filter_toAlign);
    
     // Initializing Normal Distributions Transform (NDT).
    pcl::NormalDistributionsTransform<PointNormalT, PointNormalT> ndt;

    // Setting scale dependent NDT parameters
    // Setting minimum transformation difference for termination condition.
    ndt.setEuclideanFitnessEpsilon (eps);
    // Setting maximum step size for More-Thuente line search.
    ndt.setStepSize (stepSize);
    //Setting Resolution of NDT grid structure (VoxelGridCovariance).
    ndt.setResolution (resolution);

    // Setting max number of registration iterations.
    ndt.setMaximumIterations (maxIter);

    // Setting point cloud to be aligned.
    ndt.setInputSource (filter_toAlign);
    // Setting point cloud to be aligned to.
    ndt.setInputTarget (filter_base);
    
    ndt.align (Final);
    
    std::cout << "0" << "," << eps << "," << ndt.getFitnessScore () << "," << ndt.hasConverged () << ",";
    
    return ndt.getFinalTransformation();
  }
  
  //adaptive correlation distances, not fully implemented to match:
  //Registration of point clouds using sample-sphere and adaptive distance restriction
  //DOI	 10.1007/s00371-011-0580-0
  inline Eigen::Matrix4f getICPAdaptTransform(CloudPtr base, CloudPtr toAlign)
  {
    pcl::IterativeClosestPoint<pcl::PointXYZI, pcl::PointXYZI> icp;
    pcl::PointCloud<pcl::PointXYZI> Final;
    Eigen::Matrix4f trans = Eigen::Matrix4f::Identity();
    double r_start = 1.0;
    double r_stop = 0.1;
    double r_scale = 0.75;
    double r;
    
    cout << "Warning ADAPT method not fully implemented\n";
    
    r = r_start;
    while (r > r_stop)
    {
      icp.setMaxCorrespondenceDistance (r);
      icp.setMaximumIterations (2);
      icp.setTransformationEpsilon (1e-10);
      icp.setInputSource(toAlign);
      icp.setInputTarget(base);
      icp.align(Final,trans);
      
      trans = icp.getFinalTransformation()*trans;
      r = r*r_scale;
      cout << "Step @ r=" << r <<"\n";
    }
    
    std::cout << "Converged: " << icp.hasConverged () << std::endl
		<< "Fitness: " << icp.getFitnessScore () << std::endl;
    
    return trans;
  }
  
  inline Eigen::Matrix4f getNDTTransform(CloudPtr base, CloudPtr toAlign)
  {
    // Filtering input scan to roughly 10% of original size to increase speed of registration.
    pcl::PointCloud<pcl::PointXYZI>::Ptr filtered_cloud (new pcl::PointCloud<pcl::PointXYZI>);
    
    char *stepSzC,*iter,*epsC,*resC;
    double stepSize = 0.1;
    int maxIter = 35;
    double eps = 0.01;
    double resolution = 1.0;
    
    stepSzC = getenv("STEP_SIZE");
    iter = getenv("MAX_ITER");
    epsC = getenv("CONV_EPS");
    resC = getenv("NDT_RES");
    if (resC)
    {
      resolution = atof(resC);
    }
    if (stepSzC)
    {
      stepSize = atof(stepSzC);
    }
    if (iter)
    {
      maxIter =atoi(iter);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }

    // Initializing Normal Distributions Transform (NDT).
    pcl::NormalDistributionsTransform<pcl::PointXYZI, pcl::PointXYZI> ndt;

    // Setting scale dependent NDT parameters
    // Setting minimum transformation difference for termination condition.
    ndt.setEuclideanFitnessEpsilon (eps);
    // Setting maximum step size for More-Thuente line search.
    ndt.setStepSize (stepSize);
    //Setting Resolution of NDT grid structure (VoxelGridCovariance).
    ndt.setResolution (resolution);

    // Setting max number of registration iterations.
    ndt.setMaximumIterations (maxIter);

    // Setting point cloud to be aligned.
    ndt.setInputSource (toAlign);
    // Setting point cloud to be aligned to.
    ndt.setInputTarget (base);

    // Calculating required rigid transform to align the input cloud to the target cloud.
    pcl::PointCloud<pcl::PointXYZI>::Ptr output_cloud (new pcl::PointCloud<pcl::PointXYZI>);
    ndt.align (*output_cloud);
	      
    std::cout << "0" << "," << eps << "," << ndt.getFitnessScore () << "," << ndt.hasConverged () << ",";
	      
    return ndt.getFinalTransformation();
  }

  inline Eigen::Matrix4f getPairTransform(CloudPtr base, CloudPtr toAlign)
  {
    CloudPtr src (new Cloud);
    CloudPtr tgt (new Cloud);
    pcl::ApproximateVoxelGrid<PointType> grid;
    
    double leafSize = 0.05;
    int kNorms = 10;
    char *leafStr,*kNStr,*corr,*iter,*epsC;
    double maxCorrDist = 0.5;
    int maxIter = 40;
    double eps = 1e-8;
    
    corr = getenv("CORR_DIST");
    iter = getenv("MAX_ITER");
    epsC = getenv("CONV_EPS");
    leafStr = getenv("LEAF_SIZE");
    kNStr = getenv("K_NORMS");
    if (corr)
    {
      maxCorrDist = atof(corr);
    }
    if (iter)
    {
      maxIter =atoi(iter);
    }
    if (epsC)
    {
      eps = atof(epsC);
    }
    if (leafStr)
    {
      leafSize = atof(leafStr);
    }
    if (kNStr)
    {
      kNorms = atoi(kNStr);
    }
    
    grid.setLeafSize (leafSize,leafSize,leafSize);
    grid.setInputCloud (base);
    grid.filter (*src);

    grid.setInputCloud (toAlign);
    grid.filter (*tgt);
    
    // Compute surface normals and curvature
    PointCloudWithNormals::Ptr points_with_normals_src (new PointCloudWithNormals);
    PointCloudWithNormals::Ptr points_with_normals_tgt (new PointCloudWithNormals);
    PointCloudWithNormals::Ptr reg_result (new PointCloudWithNormals);

    pcl::NormalEstimation<PointType, PointNormalT> norm_est;
    pcl::search::KdTree<PointType>::Ptr tree (new pcl::search::KdTree<PointType> ());
    norm_est.setSearchMethod (tree);
    norm_est.setKSearch (kNorms);
    
    norm_est.setInputCloud (src);
    norm_est.compute (*points_with_normals_src);
    pcl::copyPointCloud (*src, *points_with_normals_src);

    norm_est.setInputCloud (tgt);
    norm_est.compute (*points_with_normals_tgt);
    pcl::copyPointCloud (*tgt, *points_with_normals_tgt);

    //
    // Instantiate our custom point representation (defined above) ...
    MyPointRepresentation point_representation;
    // ... and weight the 'curvature' dimension so that it is balanced against x, y, and z
    float alpha[4] = {1.0, 1.0, 1.0, 1.0};
    point_representation.setRescaleValues (alpha);

    //
    // Align
    pcl::IterativeClosestPointNonLinear<PointNormalT, PointNormalT> reg;
    reg.setEuclideanFitnessEpsilon (eps);
    // Set the maximum distance between two correspondences (src<->tgt) to 10cm
    // Note: adjust this based on the size of your datasets
    reg.setMaxCorrespondenceDistance (maxCorrDist);  
    // Set the point representation
    reg.setPointRepresentation (boost::make_shared<const MyPointRepresentation> (point_representation));

    reg.setInputSource (points_with_normals_src);
    reg.setInputTarget (points_with_normals_tgt);

    // Run the same optimization in a loop and visualize the results
    reg.setMaximumIterations (maxIter);

    // Estimate
    reg.align (*reg_result);
	      
    std::cout << maxCorrDist << "," << eps << "," << reg.getFitnessScore () << "," << reg.hasConverged () << ",";

    // Get the transformation from target to source
    return reg.getFinalTransformation();
  }
};

#endif



