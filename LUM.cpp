/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2017, Jochen Sprickerhof
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

 //This function does global registration with LUM and GPS
 //the point clouds are first placed according to GPS and then registered using LUM

#include <pcl/console/parse.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/transforms.h>
#include <pcl/registration/lum.h>
#include <pcl/registration/correspondence_estimation.h>

#include <iostream>
#include <string>

#include <vector>
#include <sstream>

#include "Transformer.hpp"
#include "mapgmuutils.hpp"

//------ Data Start ------
double dist = 2.5;
int iter = 10;
int lumIter = 1;
double loopDist = 0.1;
int loopCount = 20;
std::string gps_file;
std::string fBase;
int minIdx=0, maxIdx=0;
int stepSz=1; //step size for skipping data
//------ Data End ------

//typedef pcl::PointXYZ PointType;
// typedef pcl::PointCloud<PointType> Cloud;
// typedef Cloud::ConstPtr CloudConstPtr;
// typedef Cloud::Ptr CloudPtr;
typedef std::pair<std::string, CloudPtr> CloudPair;
typedef std::vector<CloudPair> CloudVector;

void createClouds(const std::string base_name, const std::string gps_name,
	                int min, int max, int step, CloudVector& clouds)
{
	Transformer transform;
	Eigen::Matrix4f transMat = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f tmpTrans = Eigen::Matrix4f::Identity();
	Eigen::Matrix4f firstRot = Eigen::Matrix4f::Identity();
	FILE *tStampFile;
	long long tStamp;
	double currTime=0,deltaTime,startTime;

	mapgmu::PathFinder pf;

	cout<<"- Loading "<<gps_name<<endl;
	pf.loadKml(gps_name,"my path");

	char pcd_name[64];
	char out_name[64];
	char ts_name[64];

	//Load time stampfile
	for (int i=min;i<=max;i+=step)
	{
		CloudPtr cloud_ (new Cloud);
		CloudPtr transedCloud_ (new Cloud);

		//load cloud
		sprintf(pcd_name,"%s/%04d.pcd",base_name.c_str(),i);
		sprintf(out_name,"%s/out-%04d.pcd",base_name.c_str(),i);
		pcl::io::loadPCDFile(pcd_name,*cloud_);

		//load time stamp
		sprintf(ts_name,"%s/%04d.ts",base_name.c_str(),i);
		if ((tStampFile = fopen(ts_name,"r")))
		{
			fscanf(tStampFile,"%lld",&tStamp);
			if (i==min)startTime = tStamp / 1.0e6;

			deltaTime = (tStamp / 1.0e6) - (currTime + startTime);
			currTime = (tStamp / 1.0e6) - startTime;
			fclose(tStampFile);
		}
		else
		{
			cerr << "! Error: File: "<<ts_name << " not found, need timestamp\n" << endl;
			exit(1);
		}

		tmpTrans=transform.getGpsTransform(transMat,&pf,currTime-deltaTime,deltaTime);

		if(i==min)
		{
			PointXYZ p1,p2;
			p1.x=0;p1.y=0;p1.z=0;
			p2.x=0;p2.y=1;p2.z=0;
			firstRot = pf.rotateAlignment2D(p1,p2,0);
			transMat = transMat*firstRot;
		}
		else
		{
			transMat = transMat*tmpTrans;
		}

		cout << i << " " << pf.getDeltaForTime((currTime-deltaTime),deltaTime).norm() << endl;
		pcl::transformPointCloud (*cloud_,*transedCloud_,transMat);
		clouds.push_back (CloudPair (out_name, transedCloud_)); //name and point cloud
		std::cout << "- Loaded file: " << pcd_name << " size: " << transedCloud_->size () << std::endl;
	}

}

void
usage (char ** argv)
{
	std::cerr << "usage: " << argv[0] << endl
	<< "\t-fb <path-containing-pcd-files>" << endl
	<< "\t-gps <*.kml>" << endl
	<< "\t-minIdx <min-file-index>" << endl
	<< "\t-maxIdx <max-file-index>" << endl;
	return;
}

bool parse_argument(int argc, char **argv)
{
	//check must-have flags
	if ( !pcl::console::find_switch (argc, argv, "-gps") ||
	     !pcl::console::find_switch (argc, argv, "-fb") ||
			 !pcl::console::find_switch (argc, argv, "-minIdx") ||
	     !pcl::console::find_switch (argc, argv, "-maxIdx") )
	{
		usage (argv);
		return false;
	}

	pcl::console::parse_argument (argc, argv, "-d", dist);

	pcl::console::parse_argument (argc, argv, "-i", iter);

	pcl::console::parse_argument (argc, argv, "-l", lumIter);

	pcl::console::parse_argument (argc, argv, "-D", loopDist);

	pcl::console::parse_argument (argc, argv, "-c", loopCount);

	pcl::console::parse_argument (argc, argv, "-gps", gps_file);

  pcl::console::parse_argument (argc, argv, "-fb", fBase);

	pcl::console::parse_argument (argc, argv, "-minIdx", minIdx);
	pcl::console::parse_argument (argc, argv, "-maxIdx", maxIdx);


	if (pcl::console::find_switch (argc, argv, "-s") )
	{
		pcl::console::parse_argument (argc, argv, "-s", stepSz);
	}

	return true;
}


int
main (int argc, char **argv)
{

	if(parse_argument(argc,argv)==false) return 0;

  CloudVector clouds;
	createClouds(fBase, gps_file, minIdx, maxIdx, stepSz, clouds);

	//create lum
	pcl::registration::LUM<PointType> lum;
	lum.setMaxIterations (lumIter);
	lum.setConvergenceThreshold (0.001f);

	for(CloudPair& c : clouds)
	{
		lum.addPointCloud (c.second);
	}

	for (int it = 0; it < iter; it++)
	{
		for (size_t i = 1; i < clouds.size (); i++)
			for (size_t j = 0; j < i; j++)
			{
				Eigen::Vector4f ci, cj;
				pcl::compute3DCentroid (*(clouds[i].second), ci);
				pcl::compute3DCentroid (*(clouds[j].second), cj);
				Eigen::Vector4f diff = ci - cj;

				//cout<<"dist between "<<clouds[i].first<<" and "<<clouds[j].first<< " = "<<diff.norm ()<<endl;

				//connect i and j if they are adjacent in sequence or
				//they are far apart in sequence but close in space
				if( (i - j == 1) || (diff.norm () < loopDist && i - j > loopCount))
				{
					//if(i - j > loopCount)
					std::cout << "add connection between " << i << " (" << clouds[i].first << ") and " << j << " (" << clouds[j].first << ")" << std::endl;
					pcl::registration::CorrespondenceEstimation<PointType, PointType> ce;
					ce.setInputTarget (clouds[i].second);
					ce.setInputSource (clouds[j].second);
					pcl::CorrespondencesPtr corr (new pcl::Correspondences);
					ce.determineCorrespondences (*corr, dist);
					if (corr->size () > 2)
						lum.setCorrespondences (j, i, corr);
				}
			}

			lum.compute ();

			for(size_t i = 0; i < lum.getNumVertices (); i++)
			{
      //std::cout << i << ": " << lum.getTransformation (i) (0, 3) << " " << lum.getTransformation (i) (1, 3) << " " << lum.getTransformation (i) (2, 3) << std::endl;
				clouds[i].second = lum.getTransformedCloud (i);
			}
		}

		for(size_t i = 0; i < lum.getNumVertices (); i++)
		{
			std::string result_filename (clouds[i].first);
			result_filename = result_filename.substr (result_filename.rfind ("/") + 1);
			pcl::io::savePCDFileBinary (result_filename.c_str (), *(clouds[i].second));
    //std::cout << "saving result to " << result_filename << std::endl;
		}

		return 0;
	}
