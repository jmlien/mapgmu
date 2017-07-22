#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <boost/algorithm/string.hpp>
#include <boost/asio.hpp>
#include <boost/version.hpp>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <pcl/features/don.h>
#include <pcl/search/organized.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/Vertices.h>
#include "DonEstimator.hpp"

#ifndef DON_FILTER_HPP
#define DON_FILTER_HPP

using namespace std;
using namespace pcl;

//class to filter cloud based on the normal or the intensity

template<typename PointType>
class DoNFilter
{
public:
  typedef typename PointCloud<PointType>::iterator CloudItr;
  typedef PointCloud<PointType> Cloud;
  typedef typename PointCloud<PointType>::Ptr CloudPtr;
  typedef PointNormal PointNT;
  typedef PointCloud<PointNT> NCloud;
  typedef PointCloud<PointNT>::iterator NCloudItr;
  typedef NCloud::Ptr NCloudPtr;

  #define MIN_K 25
  #define MAX_K 50

  inline void filter(CloudPtr cloud_, NCloudPtr doncloud_filtered)
  {
    DonEstimator<PointType> donEst;

    // Create output cloud for DoN results
    PointCloud<PointNT>::Ptr doncloud_ (new pcl::PointCloud<PointNT>);

    // Compute DoN
    donEst.setInputCloud(cloud_);
    donEst.calculateDoN (doncloud_);
    
    pcl::ConditionAnd<PointNormal>::Ptr range_cond (
      new pcl::ConditionAnd<PointNormal> ()
      );
    range_cond->addComparison (pcl::FieldComparison<PointNormal>::ConstPtr (
				new pcl::FieldComparison<PointNormal> ("curvature", pcl::ComparisonOps::GT,0.25))
			      );
    range_cond->addComparison (pcl::FieldComparison<PointNormal>::ConstPtr (
				new pcl::FieldComparison<PointNormal> ("curvature", pcl::ComparisonOps::LT, 0.9))
			      );
    // Build the filter
    pcl::ConditionalRemoval<PointNormal> condrem (range_cond);
    condrem.setInputCloud (doncloud_);

    // Apply filter
    condrem.filter (*doncloud_filtered);
    
    //cout << "Got " << doncloud_filtered->size() <<" of " << cloud_->size() << " total Pts\n";
  }
  
  inline void intensityFilter(double numStdDev, CloudPtr cloud_in, CloudPtr cloud_out)
  {
    int ct = 0;
    double int_sum,int_sq_sum,int_mean;
    CloudItr npIt = cloud_in->begin();
    
    while (npIt<cloud_in->end())
    {
      int_sum += (*npIt).intensity;
      int_sq_sum += (*npIt).intensity*(*npIt).intensity;
      ct++;
      npIt++;
    }
    int_mean = int_sum/ct;
    int_sum = sqrt(int_sq_sum/ct - int_mean*int_mean)*numStdDev;
    
    typename pcl::ConditionAnd<PointType>::Ptr range_cond (
      new pcl::ConditionAnd<PointType> ()
      );
    
    typename pcl::FieldComparison<PointType>::ConstPtr int_comp(
	new pcl::FieldComparison<PointType> ("intensity", pcl::ComparisonOps::GT, int_mean+int_sum)
      );

    range_cond->addComparison (int_comp);
    // Build the filter
    pcl::ConditionalRemoval<PointType> condrem (range_cond);
    condrem.setInputCloud (cloud_in);

    // Apply filter
    condrem.filter (*cloud_out);
    //cout << "Got " << cloud_out->size() <<" of " << cloud_in->size() << " total Pts\n";
  }
};

#endif