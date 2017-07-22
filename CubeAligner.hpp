#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>
#include <string>
#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <Eigen/Dense>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef CUBE_ALIGNER_H
#define CUBE_ALIGNER_H

#define RAND_SCALE 100000
#define V_CT 7

#ifdef ADAPT_SIMPLEX
  #define ALPHA 1.0
  #define BETA (4.0/3.0)
  #define GAMMA (2.0/3.0)
  #define DELTA (5.0/6.0)
#else
  #define ALPHA 1.0
  #define BETA 2.0
  #define GAMMA 0.5
  #define DELTA 0.5
#endif

#define CUBE_MAX_ITER 500
//#define E_VAL 0.1
#define EPS_VAL 1200
#define MIN_CLUSTER_DENS 0.66

#ifndef MY_PI
#define MY_PI 3.141592
#endif

using namespace std;
using namespace pcl;

static bool isRandSeeded = false;

static double E_VAL = 0.2;
  
class CVertex
{
public:
  CVertex(){
    _alpha = 0;
    _beta = 0;
    _gamma = 0;
    _x = 0;
    _y = 0;
    _z = 0;
  }
  
  CVertex(double a, double b, double g, double x, double y, double z)
  {
    _alpha = a;
    _beta = b;
    _gamma = g;
    _x = x;
    _y = y;
    _z = z;
    _score = 0;
  }
  
  CVertex(Eigen::Matrix4f input)
  {
    double cb;
    _x = input(0,3);
    _y = input(1,3);
    _z = input(2,3);
    _beta = asin(-input(2,0));
    cb = cos(_beta);
    _alpha = asin(input(1,0)/cb);
    _gamma = asin(input(2,1)/cb);
    
    _score = 0;
  }
  
  /*CVertex(Eigen::Matrix4d input)
  {
    double cb;
    _x = input(0,3);
    _y = input(1,3);
    _z = input(2,3);
    _beta = asin(-input(2,0));
    cb = cos(_beta);
    _alpha = asin(input(1,0)/cb);
    _gamma = asin(input(2,1)/cb);
    
    _score = 0;
  }*/
  
#define RAD_SCALE 1.0
  
  CVertex(int i, double d, const CVertex& vert)
  {
    if (!isRandSeeded)
    {
      srand(time(NULL));
      isRandSeeded = true;
    }
    int m = i/7;
    _alpha = vert._alpha;
    _beta = vert._beta;
    _gamma = vert._gamma;
    _x = vert._x;
    _y = vert._y;
    _z = vert._z;
    
    switch(i%7)
    {
      case 0:
	_alpha += (1-2*m)*d/RAD_SCALE;
	break;
      case 1:
	_beta += (1-2*m)*d/RAD_SCALE;
	break;
      case 2:
	_gamma += (1-2*m)*d/RAD_SCALE;
	break;
      case 3:
	_x += (1-2*m)*d;
	break;
      case 4:
	_y += (1-2*m)*d;
	break;
      case 5:
	_z += (1-2*m)*d;
	break;
      default:
	break;
    }
    /*switch(j)
    {
      case 1:
	_alpha += d/10.0;
	break;
      case 2:
	_beta += d/10.0;
	break;
      case 3:
	_gamma += d/10.0;
	break;
      case 4:
	_x += d/10.0;
	break;
      case 5:
	_y += d/10.0;
	break;
      case 6:
	_z += d/10.0;
	break;
      default:
	break;
    }*/
    _score = 0;
  }
  
  CVertex(const CVertex& cg)
  {
    _alpha = cg._alpha;
    _beta = cg._beta;
    _gamma = cg._gamma;
    _x = cg._x;
    _y = cg._y;
    _z = cg._z;
    _score = 0;
  }
  
  CVertex(CVertex *cgArr, int count)
  {
    int i;
    for (i=0;i<count;i++)
    {
      _alpha = cgArr[i]._alpha/count;
      _beta = cgArr[i]._beta/count;
      _gamma = cgArr[i]._gamma/count;
      _x = cgArr[i]._x/count;
      _y = cgArr[i]._y/count;
      _z = cgArr[i]._z/count;
    }
    _score = 0;
  }
  
  CVertex(double maxRot, double maxTran)
  {
    if (!isRandSeeded)
    {
      srand(time(NULL));
      isRandSeeded = true;
    }
    _alpha = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxRot)/(double)RAND_SCALE;
    _beta = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxRot)/(double)RAND_SCALE;
    _gamma = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxRot)/(double)RAND_SCALE;
    _x = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxTran)/(double)RAND_SCALE;
    _y = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxTran)/(double)RAND_SCALE;
    _z = (((rand()% (RAND_SCALE*2))-RAND_SCALE)*maxTran)/(double)RAND_SCALE;
    _score = 0;
  }
  
  inline double getDist()
  {
    return sqrt(_x*_x + _y*_y + _z*_z);
  }
  
  inline double getMaxDiff(CVertex comp)
  {
    double diff = 0;
    
    if (fabs(_x-comp._x) > diff)
    {
      diff = fabs(_x-comp._x);
    }
    if (fabs(_y-comp._y) > diff)
    {
      diff = fabs(_y-comp._y);
    }
    if (fabs(_z-comp._z) > diff)
    {
      diff = fabs(_z-comp._z);
    }
    if (fabs(_alpha-comp._alpha) > diff)
    {
      diff = fabs(_alpha-comp._alpha);
    }
    if (fabs(_beta-comp._beta) > diff)
    {
      diff = fabs(_beta-comp._beta);
    }
    if (fabs(_gamma-comp._gamma) > diff)
    {
      diff = fabs(_gamma-comp._gamma);
    }
    
    return diff;
  }
  
  inline Eigen::Matrix4f toTranslation()
  {
    Eigen::Matrix4f outMat;
    double ca = cos(_alpha);
    double cb = cos(_beta);
    double cg = cos(_gamma);
    double sa = sin(_alpha);
    double sb = sin(_beta);
    double sg = sin(_gamma);
    
    //build the trans-rot matrix
    outMat(0,0) = ca*cb;
    outMat(0,1) = ca*sb*sg - sa*cg;
    outMat(0,2) = sa*sg - ca*sb*cg;
    outMat(1,0) = sa*cb;
    outMat(1,1) = ca*cg - sa*sb*sg;
    outMat(1,2) = sa*sb*cg - ca*sg;
    outMat(2,0) = -sb;
    outMat(2,1) = cb*sg;
    outMat(2,2) = cb*cg;
    outMat(3,0) = 0;outMat(3,1)=0;outMat(3,2)=0;
    outMat(0,3) = _x;
    outMat(1,3) = _y;
    outMat(2,3) = _z;
    outMat(3,3) = 1;
    
    return outMat;
  }
  
  inline void setScore(int score)
  {
    _score = score;
  }
  
  inline int getScore()
  {
    return _score;
  }
  
  inline static int sorter(const void *v1, const void *v2)
  {
    CVertex *cv1,*cv2;
    if (!v1 || !v2)
    {
      return 0;
    }
    cv1 = (CVertex*)v1;
    cv2 = (CVertex*)v2;
    
    return (cv1->_score < cv2->_score) - 
      2*(cv1->_score > cv2->_score);
  }
  
  inline CVertex operator+(const CVertex& vIn)
  {
    CVertex outVert(this,1);
    outVert += vIn;
    
    return outVert;
  }
  
  inline CVertex operator+=(const CVertex& vIn)
  {
    _alpha += vIn._alpha;
    _beta += vIn._beta;
    _gamma += vIn._gamma;
    _x += vIn._x;
    _y += vIn._y;
    _z += vIn._z;
    
    return *this;
  }
  
  inline CVertex operator-(const CVertex& vIn)
  {
    CVertex outVert(this,1);
    outVert -= vIn;
    
    return outVert;
  }
  
  inline CVertex operator-=(const CVertex& vIn)
  {
    _alpha -= vIn._alpha;
    _beta -= vIn._beta;
    _gamma -= vIn._gamma;
    _x -= vIn._x;
    _y -= vIn._y;
    _z -= vIn._z;
    
    return *this;
  }
  
  inline CVertex operator*(float alpha)
  {
    CVertex outVert(this,1);
    outVert *= alpha;
    
    return outVert;
  }
  
  inline CVertex operator*=(float alpha)
  {
    _alpha *= alpha;
    _beta *= alpha;
    _gamma *= alpha;
    _x *= alpha;
    _y *= alpha;
    _z *= alpha;
    
    return *this;
  }
  
  inline double getX(){return _x;}
  inline double getY(){return _y;}
  inline double getZ(){return _z;}
  
private:
  friend std::ostream& operator<<(std::ostream &strm, const CVertex& vert);
  double _alpha;
  double _beta;
  double _gamma;
  double _x;
  double _y;
  double _z;
  int _score;
};

//add the scalar*matrix operator outside the class
inline CVertex operator*(float alpha, const CVertex& vert)
{
  CVertex outMat(vert);
  outMat *= alpha;
  return outMat;
}

//String writer
std::ostream& operator<<(std::ostream &strm, const CVertex& vert)
{
  return strm << "(" << vert._x << "," << vert._y << "," << vert._z << "," << vert._alpha << "," << vert._beta << "," << vert._gamma <<")";
}


class CubeAligner
{
public:
  typedef PointCloud<PointXYZI> Cloud;
  typedef PointCloud<PointXYZI>::Ptr CloudPtr;
  typedef PointCloud<PointXYZI>::iterator CloudItr;
  
  CubeAligner(CloudPtr cloudBase, double e)
  {
    CloudItr npIt = cloudBase->begin();
    _min_x=0;
    _max_x=0;
    _min_y=0;
    _max_y=0;
    _min_z=0;
    _max_z=0;
    _arrSz = 0;
    
    while (npIt<cloudBase->end())
    {
      if (_min_x==0)
      {
	_min_x = (*npIt).x;
	_max_x = (*npIt).x;
	_min_y = (*npIt).y;
	_max_y = (*npIt).y;
	_min_z = (*npIt).z;
	_max_z = (*npIt).z;
      }
      else
      {
	if ((*npIt).x < _min_x)
	{
	  _min_x = (*npIt).x;
	}
	if ((*npIt).x > _max_x)
	{
	  _max_x = (*npIt).x;
	}
	if ((*npIt).y < _min_y)
	{
	  _min_y = (*npIt).y;
	}
	if ((*npIt).y > _max_y)
	{
	  _max_y = (*npIt).y;
	}
	if ((*npIt).z < _min_z)
	{
	  _min_z = (*npIt).z;
	}
	if ((*npIt).z > _max_z)
	{
	  _max_z = (*npIt).z;
	}
      }
      npIt++;
    }
    int gf=0;
    
    do
    {
      _e = e + ((gf++) * 0.01);
      _xCount = (_max_x - _min_x)/_e;
      _yCount = (_max_y - _min_y)/_e;
      _zCount = (_max_z - _min_z)/_e;
    }
    while (_xCount*_yCount*_zCount <= 0);
    
    //printf("Setup cube grid for (%.2f,%.2f,%.2f)->(%.2f,%.2f,%.2f)\n"
    //"%ld grid points\n",
	   //_min_x,_min_y,_min_z,_max_x,_max_y,_max_z,arraySizeNeeded());
    _vBase = (unsigned char *)malloc(arraySizeNeeded());
    memset(_vBase,0,arraySizeNeeded());
    
    
    //printf("%d/%ld filled grid points\n",cloudToArray(cloudBase,_vBase),_xCount*_yCount*_zCount);
    cloudToArray(cloudBase,_vBase);
  }
  
  ~CubeAligner()
  {
    if (_vBase)
    {
      free(_vBase);
    }
  }
  
  inline void idxToXYZ(long idx, double *x, double *y, double *z)
  {
    int ix,iy,iz;
    ix = idx%_xCount;
    iy = (idx/_xCount)%_yCount;
    iz = idx/_xCount/_yCount;
    
    *x = _min_x + ix*_e;
    *y = _min_y + iy*_e;
    *z = _min_z + iz*_e;
    
    return;
  }
    
  inline long toIdx(long ix, long iy, long iz)
  {
    if (ix < 0 || ix >= _xCount || 
      iy < 0 || iy >= _yCount || 
      iz < 0 || iz >= _zCount)
      return -1;
    
    return ix+iy*_xCount+iz*_xCount*_yCount;
  }
  
  inline size_t arraySizeNeeded()
  {
    size_t totalCount;
    
    if (!_arrSz)
    {
      totalCount = _xCount*_yCount*_zCount;
      
      totalCount += (totalCount%8)?8:0;
      
      _arrSz = totalCount/8;
    }
    
    return _arrSz;
  }
  
  inline int cloudToArray(CloudPtr cloud, unsigned char *v)
  {
    CloudItr npIt = cloud->begin();
    long ix,iy,iz;
    int bitIdx;
    long long idx;
    size_t arrIdx;
    int ptCount = 0;
    int matchCount = 0;
    
    while (npIt<cloud->end())
    {
      ptCount++;
      ix = ((*npIt).x-_min_x)/_e;
      iy = ((*npIt).y-_min_y)/_e;
      iz = ((*npIt).z-_min_z)/_e;
      //ignore points outsize range
      if (ix < 0 || iy < 0 || iz < 0
	|| ix >= _xCount || iy >= _yCount || iz >= _zCount)
      {npIt++;continue;}
      
      matchCount++;
      
      idx = ix+iy*_xCount+iz*_xCount*_yCount;
      arrIdx = idx/8;
      bitIdx = idx%8;
      v[arrIdx] |= (0x80 >> bitIdx);
      
      npIt++;
    }
    
    return matchCount;
  }
  
  class CCluster
  {
  public:
    CCluster(CubeAligner *cube){_my_cube = cube;}
    
    inline void addIdx(long idx)
    {
      _my_pts.push_back(idx);
    }
    
    inline void addIdx(long ix, long iy, long iz)
    {
      _my_pts.push_back(_my_cube->toIdx(ix,iy,iz));
    }
    
    inline void calcDimensions()
    {
      std::vector<long>::iterator it;
      double x,y,z;
      double mX=0,mY=0,mZ=0; // max values
      //set origin to obsurdly large value to fin min value
      _x_orig = 1e9;
      _y_orig = 1e9;
      _z_orig = 1e9;
      
      for (it = _my_pts.begin();it < _my_pts.end(); it++)
      {
	_my_cube->idxToXYZ(*it,&x,&y,&z);
	
	if (x < _x_orig){_x_orig = x;}
	if (y < _y_orig){_y_orig = y;}
	if (z < _z_orig){_z_orig = z;}
	if (x > mX){mX = x;}
	if (y > mY){mY = y;}
	if (z > mZ){mZ = z;}
      }
      _dx = mX - _x_orig;
      _dy = mY - _y_orig;
      _dz = mZ - _z_orig;
      
      return;
    }
    
    inline bool isInside(long idx)
    {
      double x,y,z;
      _my_cube->idxToXYZ(idx,&x,&y,&z);
      
      if (x >= _x_orig && x <= (_x_orig+_dx) &&
	y >= _y_orig && y <= (_y_orig+_dy) &&
	z >= _z_orig && z <= (_z_orig+_dz))
      {
	return true;
      }
      else
      {
	return false;
      }
    }
    
    inline bool includes(long idx)
    {
      std::vector<long>::iterator it;
      
      for (it = _my_pts.begin(); it < _my_pts.end(); it++)
      {
	if (*it == idx)
	{
	  return true;
	}
      }
      
      return false;
    }
    
    inline PointXYZI getCOG()
    {
      std::vector<long>::iterator lIt;
      PointXYZI p1;
      double x,y,z;
      p1.x = 0;
      p1.y = 0;
      p1.z = 0;
      
      for (lIt = _my_pts.begin(); lIt < _my_pts.end() ; lIt++)
      {
	_my_cube->idxToXYZ(*lIt,&x,&y,&z);
	p1.x += x/_my_pts.size();
	p1.y += y/_my_pts.size();
	p1.z += z/_my_pts.size();
      }
      
      return p1;
    }
    
    inline CloudPtr toPolygon(int i)
    {
      CloudPtr cloudOut(new Cloud);
      std::vector<long>::iterator lIt;
      
      PointXYZI p1;
      double x,y,z;
      
      for (lIt = _my_pts.begin(); lIt < _my_pts.end() ; lIt++)
      {
	_my_cube->idxToXYZ(*lIt,&x,&y,&z);
	p1.x = x;
	p1.y = y;
	p1.z = z;
  p1.intensity = i;
	cloudOut->push_back(p1);
      }
      
//       PointXYZI p1,p2,p3,p4,p5,p6,p7,p8;
//       p1.x = _x_orig;
//       p1.y = _y_orig;
//       p1.z = _z_orig;
//       p2.x = _x_orig+_dx;
//       p2.y = _y_orig;
//       p2.z = _z_orig;
//       p3.x = _x_orig+_dx;
//       p3.y = _y_orig+_dy;
//       p3.z = _z_orig;
//       p4.x = _x_orig;
//       p4.y = _y_orig+_dy;
//       p4.z = _z_orig;
//       p5.x = _x_orig;
//       p5.y = _y_orig;
//       p5.z = _z_orig + _dz;
//       p6.x = _x_orig+_dx;
//       p6.y = _y_orig;
//       p6.z = _z_orig + _dz;
//       p7.x = _x_orig+_dx;
//       p7.y = _y_orig+_dy;
//       p7.z = _z_orig + _dz;
//       p8.x = _x_orig;
//       p8.y = _y_orig+_dy;
//       p8.z = _z_orig + _dz;
//       cloudOut->push_back(p1);
//       cloudOut->push_back(p2);
//       cloudOut->push_back(p3);
//       cloudOut->push_back(p4);
//       cloudOut->push_back(p5);
//       cloudOut->push_back(p6);
//       cloudOut->push_back(p7);
//       cloudOut->push_back(p8);
      
      return cloudOut;
    }
    
    inline int getSize(){return _my_pts.size();}
    
    inline double getDensity()
    {
      return ((double)_my_pts.size())/((_dx*_dy*_dz)/(_my_cube->getE3()));
    }
    
    inline CCluster splitMe()
    {
      CCluster newClust(_my_cube);
      CCluster backFillClust(_my_cube);
      
      std::vector<long>::iterator cIt;
      long ix,iy,iz;
      long z;
      long hitCount,lastHitCount = -1;
      double lastDens = 0;
      
      cout << "Trying to split" << endl;
      
      for (iz=0;iz<_my_cube->getZCount();iz++)
      {
	/*hitCount = 0;
	for (ix=0;ix<_my_cube->getXCount();ix++)
	{
	  for (iy=0;iy<_my_cube->getYCount();iy++)
	  {
	    if (_my_cube->isFilledIdx(_my_cube->toIdx(ix,iy,iz)))
	    {
	      hitCount++;
	    }
	  }
	}*/
	
	for (cIt = _my_pts.begin();cIt < _my_pts.end();
	      cIt++)
	{
	  z = (*cIt)/_my_cube->getXCount()/_my_cube->getYCount();
	  if (z == iz)
	  {
	    newClust.addIdx(*cIt);
	  }
	}
	if (iz > 1)
	{
	  newClust.calcDimensions();
	  if (lastDens > 0 && newClust.getDensity() < lastDens)
	  {
	    break;
	  }
	  lastDens = newClust.getDensity();
	}
      }
      
      newClust._my_pts.clear();
      if (iz < _my_cube->getZCount())
      {
	for (cIt = _my_pts.begin();cIt < _my_pts.end();
	      cIt++)
	{
	  z = (*cIt)/_my_cube->getXCount()/_my_cube->getYCount();
	  if (z < iz)
	  {
	    newClust.addIdx(*cIt);
	  }
	  else
	  {
	    backFillClust.addIdx(*cIt);
	  }
	}
	_my_pts.clear();
	for (cIt = backFillClust._my_pts.begin();cIt < backFillClust._my_pts.end();
	      cIt++)
	{
	  addIdx(*cIt);
	}
	backFillClust._my_pts.clear();
      }
      
      return newClust;
    }
    
    inline double getXMin(){return _x_orig;}
    inline double getYMin(){return _y_orig;}
    inline double getZMin(){return _z_orig;}
    inline double getXMax(){return _x_orig+_dx;}
    inline double getYMax(){return _y_orig+_dy;}
    inline double getZMax(){return _z_orig+_dz;}
    inline double getDX(){return _dx;}
    inline double getDY(){return _dy;}
    inline double getDZ(){return _dz;}
    
    inline std::vector<long> *getIdxs(){return &_my_pts;}
    
  private:
    std::vector<long> _my_pts;
    CubeAligner *_my_cube;
    double _x_orig;
    double _y_orig;
    double _z_orig;
    double _dx;
    double _dy;
    double _dz;
  };
  
  inline double getE3(){return _e*_e*_e;}
  
  inline void removeCluster(CCluster *cluster)
  {
    std::vector<long>::iterator clustIter;
    long idx,arrIdx;
    int bitIdx;
    
    for (clustIter = cluster->getIdxs()->begin();
	 clustIter < cluster->getIdxs()->end();
	 clustIter++)
    {
      idx = *clustIter;
    
      arrIdx = idx/8;
      bitIdx = idx%8;
      _vBase[arrIdx] ^= (0x80 >> bitIdx);
    }
    
    return;
  }
  
  inline void identifyClusters()
  {
    long ix,iy,iz,idx;
    std::vector<CCluster>::iterator clustIter;
    
    _clusters.clear();
    
    for (iz=0;iz<(-_min_z)/_e;iz++)
    {
      for (ix=0;ix<_xCount;ix++)
      {
	for (iy=0;iy<_yCount;iy++)
	{
	  idx = toIdx(ix,iy,iz);
	  if ( isFilledIdx(idx) )
	  {
	    bool found = false;
	    for (clustIter = _clusters.begin(); clustIter < _clusters.end();
		 clustIter++)
	    {
	      if ((*clustIter).includes(idx))
	      {
		found = true;
		break;
	      }
	    }
	    if (found)
	    {
	      continue;
	    }
	    CCluster newClust((CubeAligner*)this);
	    //growClusterRadial(idx,idx,&newClust);
	    growCluster(idx,&newClust);
	    if (newClust.getSize() > 25)
	    {
	      newClust.calcDimensions();
	      _clusters.push_back(newClust);
	    }
	  }
	}
      }
    }
    
    splitClusters();
    
    return;
    
    for (clustIter = _clusters.begin(); clustIter < _clusters.end();
	  clustIter++)
    {
      //if ((*clustIter).getZMin() < 0)
      //{
	removeCluster(&(*clustIter));
      //}
    }
    
    _clusters.clear();
    
    for (ix=0;ix<_xCount;ix++)
    {
      for (iy=0;iy<_yCount;iy++)
      {
	for (iz=0;iz<_zCount;iz++)
	{
	  idx = toIdx(ix,iy,iz);
	  if ( isFilledIdx(idx) )
	  {
	    bool found = false;
	    for (clustIter = _clusters.begin(); clustIter < _clusters.end();
		 clustIter++)
	    {
	      if ((*clustIter).includes(idx))
	      {
		found = true;
		break;
	      }
	    }
	    if (found)
	    {
	      continue;
	    }
	    CCluster newClust((CubeAligner*)this);
	    growCluster(idx,&newClust);
	    if (newClust.getSize() > 20)
	    {
	      newClust.calcDimensions();
	      //cout << "Added Cluster -> size = " << newClust.getSize() << endl;
	      found = false;
	      for (clustIter = _clusters.begin(); clustIter < _clusters.end();
		  clustIter++)
	      {
		if (
		  (*clustIter).getDensity() < newClust.getDensity() )
		{
		  _clusters.insert(clustIter,newClust);
		  found = true;
		  break;
		}
	      }
	      if (!found)
	      {
		_clusters.push_back(newClust);
	      }
	    }
	  }
	}
      }
    }
    
    cout << "Found " << _clusters.size() << " Clusters" << endl;
    
    return;
  }
  
  inline std::vector<CCluster> *getClusters()
  {
    return &_clusters;
  }
  
  inline void splitClusters()
  {
    std::vector<CCluster>::iterator clustIter;
    for (clustIter = _clusters.begin(); clustIter < _clusters.end();
	clustIter++)
    {
      if ((*clustIter).getDensity() < MIN_CLUSTER_DENS )
      {
	_clusters.push_back((*clustIter).splitMe());
      }
    }
    
    return;
  }
  
  inline void growClusterRadial(long startIdx, long thisIdx, CCluster *clust)
  {
    long ix,iy,iz,ix2,iy2,iz2;
    int dx,dy,dz;
    long newIdx;
    bool isThisFilled;
    
    ix = startIdx%_xCount;
    iy = (startIdx/_xCount)%_yCount;
    iz = startIdx/_xCount/_yCount;
    
    isThisFilled = isFilledIdx(thisIdx);
    if (isThisFilled)
    {
      if (!clust->includes(thisIdx))
      {
	clust->addIdx(thisIdx);
      }
      else
      {
	return;
      }
    }
    
    if (startIdx==thisIdx)
    {
      for (dx=-1;dx<=1;dx++)
      {
	for (dy=-1;dy<=1;dy++)
	{
	  if (!dx && !dy){continue;}
	  newIdx = toIdx(ix+dx,iy+dy,iz);
	  if (newIdx >=  0)
	    growClusterRadial(startIdx,newIdx,clust);
	}
      }
    }
    else
    {
      ix2 = thisIdx%_xCount;
      iy2 = (thisIdx/_xCount)%_yCount;
      iz2 = thisIdx/_xCount/_yCount;
      dx = ix2-ix;
      dy = iy2-iy;
      
      //cout << "checking " << ix2<<","<<iy2<<","<<iz2<<endl;
      //cout << "ctr = " << ix<<","<<iy<<","<<iz<<endl;
      
      if (dx)
      {
	newIdx = toIdx(ix2+(dx<0?-1:1),iy2,iz2);
	if ( newIdx >= 0 && (isFilledIdx(newIdx) || isThisFilled))
	  growClusterRadial(startIdx,newIdx,clust);
	/*if (!isFilledIdx(newIdx) && isThisFilled)
	{
	  newIdx = toIdx(ix2+(dx<0?-1:1),iy2,iz2+1);
	  if (isFilledIdx(newIdx) && !clust->includes(newIdx))
	    growClusterRadial(startIdx,newIdx,clust);
	  newIdx = toIdx(ix2+(dx<0?-1:1),iy2,iz2-1);
	  if (isFilledIdx(newIdx) && !clust->includes(newIdx))
	    growClusterRadial(startIdx,newIdx,clust);
	}*/
      }
      if (dy)
      {
	newIdx = toIdx(ix2,iy2+(dy<0?-1:1),iz2);
	if ( newIdx >= 0 && (isFilledIdx(newIdx) || isThisFilled))
	  growClusterRadial(startIdx,newIdx,clust);
	/*if (!isFilledIdx(newIdx) & isThisFilled)
	{
	  newIdx = toIdx(ix2,iy2+(dy<0?-1:1),iz2+1);
	  if (isFilledIdx(newIdx) && !clust->includes(newIdx))
	    growClusterRadial(startIdx,newIdx,clust);
	  newIdx = toIdx(ix2,iy2+(dy<0?-1:1),iz2-1);
	  if (isFilledIdx(newIdx) && !clust->includes(newIdx))
	    growClusterRadial(startIdx,newIdx,clust);
	}*/
      }
      if (dx && dy)
      {
	newIdx = toIdx(ix2+(dx<0?-1:1),iy2+(dy<0?-1:1),iz2);
	if ( newIdx >= 0 && (isFilledIdx(newIdx) || isThisFilled))
	  growClusterRadial(startIdx,newIdx,clust);
      }
    }
    
    return;
  }
  
  inline void growCluster(long startIdx, CCluster *clust)
  {
    long ix,iy,iz,idx;
    ix = startIdx%_xCount;
    iy = (startIdx/_xCount)%_yCount;
    iz = startIdx/_xCount/_yCount;
    int offset;
    
    clust->addIdx(startIdx);
    
    for (offset = -1; offset <= 1; offset += 2)
    {
      idx = toIdx(ix,iy+offset,iz);
      if (idx >= 0 && isFilledIdx(idx) && !clust->includes(idx))
      {
	growCluster(idx,clust);
      }
      else if ((idx = toIdx(ix,iy+offset*2,iz)) >= 0 && isFilledIdx(idx) &&
	!clust->includes(idx))
      {
	growCluster(idx,clust);
      }
      
      idx = toIdx(ix,iy,iz+offset);
      if (idx >= 0 && isFilledIdx(idx) && !clust->includes(idx))
      {
	growCluster(idx,clust);
      }
      else if ( (idx = toIdx(ix,iy,iz+offset*2)) >= 0 && isFilledIdx(idx) &&
      !clust->includes(idx))
      {
      growCluster(idx,clust);
      }
      
      //if (abs(ixStart - (ix+offset)) <= 3)
      //{
      idx = toIdx(ix+offset,iy,iz);
      if (idx >= 0 && isFilledIdx(idx) && !clust->includes(idx))
      {
	growCluster(idx,clust);
      }
      else if ((idx = toIdx(ix+offset*2,iy,iz)) >= 0 && isFilledIdx(idx) &&
	!clust->includes(idx))
      {
	growCluster(idx,clust);
      }
      //}
    }
    
    return;
  }
  
  inline void filterCloudCubic(CloudPtr cloud)
  {
    CloudItr npIt = cloud->begin();
    CloudPtr cloud2(new Cloud);
    ExtractIndices<PointXYZI> iextract(false);
    boost::shared_ptr<vector<int> > removeIndices(new vector<int>);
    long ix,iy,iz;
    int bitIdx;
    long long idx;
    size_t arrIdx;
    int i=0;
    int dropCount=0;
    
    removeIndices->reserve(cloud->size());
    
    while (npIt<cloud->end())
    {
      ix = ((*npIt).x-_min_x)/_e;
      iy = ((*npIt).y-_min_y)/_e;
      iz = ((*npIt).z-_min_z)/_e;
      
      if (isFilledIdx(toIdx(ix,iy,iz)))
      {
	removeIndices->push_back(i);
	dropCount++;
      }
      
      npIt++;
      i++;
    }
    //cout << "before " << cloud->size();
    iextract.setInputCloud(cloud);
    iextract.setNegative(true);
    iextract.setIndices(removeIndices);
    iextract.filter(*cloud2);
    cloud->swap(*cloud2);
    //cout << "after " << cloud->size() << endl;

    //cout << "Dropped " << dropCount << " non cubic points" << endl;
    return;
  }
  
  inline Eigen::Matrix4f alignToMatrix(Eigen::Matrix4f mat, CloudPtr cloud,
    double *eVal)
  {
    CloudPtr transCloud;
    CVertex x0;
    CVertex xr;
    CVertex xe;
    CVertex vertArr[V_CT];
    int i,j;
    int lCount=0;
    int cubeCount;
    int loopCount = 0;
    unsigned char *vArr;
    vArr = (unsigned char *)malloc(arraySizeNeeded());
    //try an initial offset in the direction of carts motion...
    CVertex v0(mat);
    double newScore;
    int didReduct;
    
    
    while(1)
    {
      loopCount++;
      for (i=0;i<V_CT;i++)
      {
	vertArr[i] = CVertex(i,*eVal,v0);
	transCloud = CloudPtr(new Cloud);
	transformPointCloud(*cloud,*transCloud,vertArr[i].toTranslation());
	memset(vArr,0,arraySizeNeeded());
	if (i==0){
	  cubeCount = cloudToArray(transCloud,vArr);
	}
	else
	{
	  cloudToArray(transCloud,vArr);
	}
  
	vertArr[i].setScore(getHitCount(vArr));
      }
      qsort(vertArr,V_CT,sizeof(CVertex),CVertex::sorter);
      didReduct = 0;
      
      while (1)
      {
  
	/*for (i=0;i<V_CT;i++)
	{
	  printf("%d,",vertArr[i].getScore());
	}
	printf("\n");*/
	for (i=0;i<V_CT;i++)
	{
	  //printf("%d,",vertArr[i].getScore());
	}
	//printf("\n");
	
	//calculate center of mass and reflection point
	x0 = CVertex(vertArr,V_CT-1);
	xr = x0 + ALPHA*(x0-vertArr[V_CT-1]);
	
	//calculate reflection points score
	transCloud = CloudPtr(new Cloud);
	transformPointCloud(*cloud,*transCloud,xr.toTranslation());
	memset(vArr,0,arraySizeNeeded());
	cloudToArray(transCloud,vArr);
	xr.setScore(getHitCount(vArr));
	
	for (i=0;i<V_CT;i++)
	{
	  if (xr.getScore() > vertArr[i].getScore())
	  {
	    break;
	  }
	}
	
	//replace the worst point with reflected and continue
	if ((i > 0 && i < V_CT-1) || xr.getScore() == vertArr[0].getScore())
	{
	  vertArr[V_CT-1] = xr;
	}
	else if (i==0)
	{
	  xe = x0 + BETA*(xr-x0);
	  //calculate expansion points score
	  transCloud = CloudPtr(new Cloud);
	  transformPointCloud(*cloud,*transCloud,xe.toTranslation());
	  memset(vArr,0,arraySizeNeeded());
	  cloudToArray(transCloud,vArr);
	  xe.setScore(getHitCount(vArr));
	  if (xe.getScore() > xr.getScore())
	  {
	    vertArr[V_CT-1] = xe;
	  }
	  else
	  {
	    vertArr[V_CT-1] = xr;
	  }
	}
	else if (i == (V_CT - 1))
	{
	  //outside contraction
	  xe = x0 + GAMMA*(vertArr[V_CT-1]-x0);
	  transCloud = CloudPtr(new Cloud);
	  transformPointCloud(*cloud,*transCloud,xe.toTranslation());
	  memset(vArr,0,arraySizeNeeded());
	  cloudToArray(transCloud,vArr);
	  xe.setScore(getHitCount(vArr));
	  if (xe.getScore()  > xr.getScore())
	  {
	    vertArr[V_CT-1] = xe;
	  }
	  else
	  {
	    i=-1;
	  }
	}
	else
	{
	  //inside contraction
	  xe = x0 - GAMMA*(vertArr[V_CT-1]-x0);
	  transCloud = CloudPtr(new Cloud);
	  transformPointCloud(*cloud,*transCloud,xe.toTranslation());
	  memset(vArr,0,arraySizeNeeded());
	  cloudToArray(transCloud,vArr);
	  xe.setScore(getHitCount(vArr));
	  if (xe.getScore()  > xr.getScore())
	  {
	    vertArr[V_CT-1] = xe;
	  }
	  else
	  {
	    i=-1;
	  }
	}
	
	if (i == -1)
	{
	  //reduction
	  if (didReduct == 5)
	  {
	    lCount = CUBE_MAX_ITER;
	  }
	  else
	  {
	    //reduce
	    for (i=1;i<V_CT;i++)
	    {
	      vertArr[i] = vertArr[0] + DELTA*(vertArr[i]-vertArr[0]);
	      transCloud = CloudPtr(new Cloud);
	      transformPointCloud(*cloud,*transCloud,vertArr[i].toTranslation());
	      memset(vArr,0,arraySizeNeeded());
	      cloudToArray(transCloud,vArr);
	      vertArr[i].setScore(getHitCount(vArr));
	    }

	    //cout << "REDUCED...";
	    didReduct++;
	  }
	}
	
	newScore = vertArr[V_CT-1].getScore();
	//printf("Next score = %d\n",vertArr[V_CT-1].getScore());
	
	qsort(vertArr,V_CT,sizeof(CVertex),CVertex::sorter);
	if (didReduct && i == V_CT)
	{
	  //dont quit early from std dev right after a reduction
	  newScore = vertArr[0].getScore();
	}
	
	//terminate the loop
	if (++lCount > CUBE_MAX_ITER)
	{
	  break;
	}
	if (newScore != vertArr[0].getScore() && 
	  getStdDev(vertArr,V_CT) < (vertArr[0].getScore()*.01))
	{
	  break;
	}
      }
      
      /*cout << "TESTING 123 " << loopCount << endl << endl << endl;
      
      //check the distance to the new and restart maybe??
      if (loopCount <= 3 && vertArr[0].getScore() < cubeCount*0.2)
      {
	lCount = 0;
	v0 = vertArr[0];
      }
      else
      {*/
	break;
      //}
    }
    
    x0 = vertArr[0];
    
    /*newScore = x0.getScore();
    for (j=0;j<10;j++)
    {
      lCount = 0;
      for (i=0;i<V_CT;i++)
      {
	xe = CVertex(i,E_VAL/100.0,x0);
	transCloud = CloudPtr(new Cloud);
	transformPointCloud(*cloud,*transCloud,xe.toTranslation());
	memset(vArr,0,arraySizeNeeded());
	cloudToArray(transCloud,vArr);
	xe.setScore(getHitCount(vArr));
	if (xe.getScore() > newScore)
	{
	  lCount = i;
	  newScore = xe.getScore();
	  xr = xe;
	}
	cout << "Offset " << i << " = " << xe.getScore() << endl;
      }
      if (lCount == 0)
      {
	break;
      }
      else
      {
	x0=xr;
      }
    }*/
    
    //E_VAL = ( E_VAL + x0.getMaxDiff(v0)*2 )/2.0;
    //if (E_VAL<0.05)E_VAL = 0.05;
    
    //*eVal = (*eVal + v0.getMaxDiff(x0))*0.5;
    //if (*eVal<0.01) *eVal = 0.01;
   
    //cout << "New E_VAL = " << *eVal << endl;
    //cout << "Starting point found:" << endl << v0 << endl;
    //cout << "Best trans = " << x0 << endl;
    
#ifdef ADAPT_SIMPLEX
    //cout << "ADAPTING!!>>!!!!!" << endl;
#endif
    
    free(vArr);
    
    return x0.toTranslation();
  }
  
  inline double getNextEval(Eigen::Matrix4f *guess, Eigen::Matrix4f *align)
  {
    CVertex c1(*guess);
    CVertex c2(*align);
    
    return c1.getMaxDiff(c2);
  }
  
  inline int getXCount(){return _xCount;}
  inline int getYCount(){return _yCount;}
  inline int getZCount(){return _zCount;}
  
private:
  
  inline CVertex findStartingPoint(CloudPtr cloud, double minX,
    double maxX, double minY, double maxY, double minZ, double maxZ,
    double stepSize)
  {
    CVertex best;
    CVertex test;
    CloudPtr transCloud;
    double x,y,z;
    unsigned char *vArr;
    vArr = (unsigned char *)malloc(arraySizeNeeded());
    
    best.setScore(0);
    
    for (x = minX;x<=maxX;x+=stepSize)
    {
      for (y = minY;y<=maxY;y+=stepSize)
      {
	for (z = minZ; z <= maxZ; z+=stepSize)
	{
	  test = CVertex(0,0,0,x,y,z);
	  cout << "testing " << test <<"..." << endl;
	  transCloud = CloudPtr(new Cloud);
	  transformPointCloud(*cloud,*transCloud,test.toTranslation());
	  memset(vArr,0,arraySizeNeeded());
	  cloudToArray(transCloud,vArr);
	  test.setScore(getHitCount(vArr));
	  if (test.getScore() > best.getScore())
	  {
	    best = test;
	  }
	}
      }
    }
    
    free(vArr);
    return best;
  }
  
  inline int getHitCount(unsigned char *v2)
  {
    int count = 0;
    int i,j;
    int byteCount = arraySizeNeeded();
    unsigned char cpByte;
    
    for (i=0;i<byteCount;i++)
    {
      cpByte = _vBase[i]&v2[i];
      for (j=0;j<8;j++)
      {
	count += ((cpByte >> j) & 0x01);
      }
    }
    
    return count;
  }
  
  inline bool isFilledIdx(long idx)
  {
    long arrIdx;
    int bitIdx;
    
    //catch out of range from toIdx()
    if (idx==-1)
    {
      return false;
    }
    
    arrIdx = idx/8;
    bitIdx = idx%8;
    return (_vBase[arrIdx] & (0x80 >> bitIdx) );
  }
    
  
  inline double getStdDev(CVertex *vertArr, int count)
  {
    double mean=0;
    double stdDev=0;
    int i;
    for (i=0;i<count;i++)
    {
      mean += vertArr[i].getScore();
    }
    mean /= count;
    
    for (i=0;i<count;i++)
    {
      stdDev += pow(vertArr[i].getScore()-mean,2.0)/(count-1);
    }
    
    stdDev = sqrt(stdDev);
    //printf("Std dev = %.2f\n",stdDev);
    return stdDev;
  }
  
  std::vector<CCluster> _clusters;
  double _min_x;
  double _max_x;
  double _min_y;
  double _max_y;
  double _min_z;
  double _max_z;
  double _e;
  long _xCount;
  long _yCount;
  long _zCount;
  size_t _arrSz;
  unsigned char* _vBase;
};




#endif
