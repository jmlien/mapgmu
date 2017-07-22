#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/Vertices.h>
#include <vector>
#include <string>
#include <cmath>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>

//cheap attempt at meshing a cloud, used in SimplePCDViewer

template<typename PointType>
class MyMeshBuilder
{
  
public:
  MyMeshBuilder(double mDist):
    tree(new TreeType),
    meshTree(new TreeType),
    minDist(mDist),
    maxDot(0.88),
    shortLen(2.0),
    dScale(2)
  {
    
  }
  typedef typename pcl::PointCloud<PointType>::Ptr CloudPtr;
  typedef typename pcl::KdTreeFLANN<PointType> TreeType;
  typedef typename TreeType::Ptr TreePtr;

  inline std::vector<pcl::Vertices> buildMesh(CloudPtr baseCloud)

  {
    std::vector<pcl::Vertices> meshVerts(0);
    std::vector<int> indices(10);
    std::vector<float> dists(10);
    float ab,adotb,x0,y0,z0;
    int nPts,i,idc2;
    PointType p2, p3;
    
    meshTree->setInputCloud(baseCloud);
    
    typename pcl::PointCloud<PointType>::iterator pIt = baseCloud->begin();
    while (pIt<baseCloud->end())
    {
      nPts = meshTree->radiusSearch(*pIt,dScale*minDist,indices,dists,10);
      if (nPts > 3)
      {
	x0 = (*pIt).x; y0 = (*pIt).y; z0 = (*pIt).z;
	p2 = baseCloud->at(indices[1]);
	for (i=2;i<nPts;i++)
	{
	  ab = dists[1]*dists[i];
	  p3 = baseCloud->at(indices[i]);
	  adotb = (p2.x - x0)*(p3.x - x0)+
		  (p2.y - y0)*(p3.y - y0)+
		  (p2.z - z0)*(p3.z - z0);
	  if (fabs(adotb/ab) < maxDot )
	  {
	    i++;
	    break;
	  }
	}
	idc2 = i-1;
      }
      else if (nPts == 3)
      {
	idc2 = 2;
      }
      else
      {
	pIt++;
	continue;
      }
      pcl::Vertices verts;
      verts.vertices.push_back( indices[0] );
      verts.vertices.push_back( indices[1] );
      verts.vertices.push_back( indices[idc2] );
      meshVerts.push_back(verts);
      
      if (idc2<nPts-1)
      {
	//How do I check if theres more??
      }
      pIt++;
    }
    
    return meshVerts;
  }
  
  TreePtr tree;
  TreePtr meshTree;
  float minDist;
  float maxDot;
  float shortLen;
  int dScale;

};