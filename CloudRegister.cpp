#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/console/parse.h>
#include <string>
#include <pcl/io/boost.h>
#include <pcl/common/transforms.h>
#include <pcl/io/pcd_io.h>
#include "Transformer.hpp"
#include <stdio.h>
#include "mapgmuutils.hpp"

#include <time.h>

//supports up to MAX_CLOCKS concurrent timers
#define MAX_CLOCKS 10
clock_t c_starts[MAX_CLOCKS];
double tmpClock[MAX_CLOCKS];

#define START_TIMER(X) do{c_starts[X] = clock();}while(0)
#define STOP_TIMER(X,CHAR_DEST) do{ \
	tmpClock[X]  = (double) (clock() - c_starts[X])/((double)CLOCKS_PER_SEC); \
	sprintf(CHAR_DEST,"%f",tmpClock[X]);}while(0)

	using namespace std;
	using namespace pcl;
	using namespace pcl::console;

	typedef pcl::PointXYZI PointType;
	typedef PointCloud<PointType> Cloud;
	typedef PointCloud<PointType>::Ptr CloudPtr;

	FILE *transFile;
	std::string transformFileName;


	void filterIntensity(int minInt,CloudPtr cloud2_)
	{
		CloudPtr cloudSwap_(new Cloud);

		pcl::ConditionAnd<PointType>::Ptr range_cond (
			new pcl::ConditionAnd<PointType> ()
		);
		range_cond->addComparison (pcl::FieldComparison<PointType>::ConstPtr (
			new pcl::FieldComparison<PointType> ("intensity", pcl::ComparisonOps::LT, minInt))
		);
		// Build the filter
		pcl::ConditionalRemoval<PointType> condrem (range_cond);
		condrem.setInputCloud (cloud2_);

		// Apply filter
		condrem.filter (*cloudSwap_);
		cout << "Intensity Filter removed " << (cloud2_->size() - cloudSwap_->size()) << " Points\n";
		cloud2_->swap(*cloudSwap_);
	}

	void SaveTransform(char *fName, Eigen::Matrix4f trans)
	{
		int i,j;
		if (!transFile)
		{
			transFile = fopen(transformFileName.c_str(),"w");
		}

		fprintf(transFile,"%s\n",fName);
		for (i=0;i<4;i++)
		{
			for (j=0;j<4;j++)
			{
				fprintf(transFile,"%.16e ",trans(i,j));
			}
			fprintf(transFile,"\n");
		}

		fflush(transFile);
	}

	void SaveFinalTransform(const char *fName, Eigen::Matrix4f trans)
	{
		int i,j;
		FILE *finalFile;

		finalFile = fopen(fName,"w");

		for (i=0;i<4;i++)
		{
			for (j=0;j<4;j++)
			{
				fprintf(finalFile,"%.16e ",trans(i,j));
			}
			fprintf(finalFile,"\n");
		}

		fclose(finalFile);
	}

	Eigen::Matrix4f LoadFinalTransform(const char *fName)
	{
		int i,j;
		FILE *finalFile;
		Eigen::Matrix4f tempMat;
		double tmpDbl;

		finalFile = fopen(fName,"r");
		if (!finalFile)
		return Eigen::Matrix4f::Identity();

		for (i=0;i<4;i++)
		{
			for (j=0;j<4;j++)
			{
				fscanf(finalFile,"%le ",&tmpDbl);
				tempMat(i,j) = tmpDbl;
			}
		}


		fclose(finalFile);

		return tempMat;
	}

	std::string transTypeStrings[] = {
		"NONE",
		"ICP",
		"NDT",
		"ICP_NL_NORM",
		"ADAPT_ICP",
		"DoN",
		"NDTDoN",
		"INT",
		"VOX",
		"CUBE",
		"C-ICP",
		"C-GPS-ICP",
		"C-GPS",
		"GPS-ICP",
		"GPS-ICP-C",
		"GPS-ICP-NC",
		"GPS",
		"C-GPS-ICP-ADAPT"
	};

	std::string transDescStrings[] = {
		"Identity Transformation",
		"Iterative Closest Point",
		"Normal Distributions Transform",
		"ICP Non Linear with Normals",
		"Adaptive Correspondance distance ICP",
		"Density of Normal filtered ICP (See DoNFilter.hpp)",
		"Density of Normal filteree DON (See DoNFilter.hpp)",
		"Intensity filtered ICP followed by full ICP",
		"Voxel filtered ICP",
		"Coarse Binary Cube Alignment with cubic filtered alignement cloud",
		"CUBE (filtered) followed by ICP",
		"CUBE (filtered) seeded by gps velocity vector followed by cubic filtered ICP then full cloud ICP",
		"CUBE (unfiltered) seeded by gps velocity vector",
		"GPS seeded ICP on cubic filtered alignment cloud then full cloud ICP",
		"GPS-ICP followed by CUBE (non-filtered)",
		"GPS seeded ICP, two passes on full clouds",
		"Tranlation only GPS velocity vector transforms",
		"C-GPS-ICP with adaptive cube size based on gps velocity estimates"
	};

	int getTranstype(std::string transStr)
	{
		int transType = 0;
		for (int i=0;i<sizeof(transTypeStrings)/sizeof(std::string);i++)
		{
			if (!transStr.compare(transTypeStrings[i]))
			{
				transType = i;
				break;
			}
		}
		return transType;
	}

	void print_trans_list()
	{
		for (int i=0;i<sizeof(transTypeStrings)/sizeof(std::string);i++)
		{
			cout << "\t" << transTypeStrings[i] << ": " << transDescStrings[i] << endl;
		}
		cout << "See Transformer.hpp for Environment Variables to overide default configuration values" << endl;
	}

	void
	usage (char ** argv)
	{
		cout << "usage: " << argv[0] << endl
		<< "\t-fb <path-containing-pcd-files>" << endl
		<< "\t-minIdx <min-file-index>" << endl
		<< "\t-maxIdx <max-file-index>" << endl
		<< "\t-t <transform_type from: ";

		for (int i=0;i<sizeof(transTypeStrings)/sizeof(std::string);i++)
		{
			cout << transTypeStrings[i] << ",";
		}
		cout << ">" << endl
		<< "\t[-trFile <output-transform-file>(defaut=trans.mat)]" << endl
		<< "\t[-s <step-count(eg 2 registers every other cloud)>]" << endl
		<< "\t[-gps <path_to_gps_kml_file, required for GPS transforms>]" << endl
		<< "\t[-sft <final-transform-output-file>] (Saves Final Transform to be loaded with lft)" << endl
		<< "\t[-lft <first-transform-input-file>] (Load First Tranform from sft in previous run)" << endl
		<< "Files are expected to have the following format {fb}/####.pcd (zero padded)"
		<< endl;
		cout << argv[0] << " -h | --help : shows this help" << endl;
		cout << argv[0] << " -trans_list : shows tranform descriptions" << endl;
		return;
	}

	int main (int argc, char ** argv)
	{
		CloudPtr cloud_ (new Cloud);
		CloudPtr cloud2_(new Cloud);
		CloudPtr finalCloud(new Cloud);
		Cloud aligned;
		int transType = 0;
		char timerString[50];
		char longTimerString[50];

		std::string fBase, minIdx, maxIdx, transStr,stepSz,intStr;
		int min,max, i,step=1;
		char fName[50];
		char fNameB[50];

		Eigen::Matrix4f transMat = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f tmpTrans = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f tmpTransB = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f firstRot = Eigen::Matrix4f::Identity();
		bool hasOutfile = false;

		CloudPtr saveCloud(new Cloud);
		CloudPtr copyCloud(new Cloud);

		mapgmu::PathFinder pf;
		bool hasGps=false;
		std::string gpsKmlFile, ftFile;
		int segIdx = 0;
		PointXYZ currPoint,prevPoint;
		CVertex currVert,prevVert,prevVel,currVel,currAcc;
		currPoint.x=0;currPoint.y=0;currPoint.z=0;
		double cDist = 0.15;
		double cSize = 0.15;
		FILE *tStampFile;
		long long tStamp;
		double currTime,deltaTime,startTime;

		int minInt;
		bool iFilter = false;

		Transformer transform;

		if (find_switch (argc, argv, "-trans_list"))
		{
			print_trans_list();
			return(0);
		}

		if (find_switch (argc, argv, "-h") ||
		find_switch (argc, argv, "--help") ||
		!( find_switch (argc, argv, "-fb") &&
		find_switch (argc, argv, "-minIdx") &&
		find_switch (argc, argv, "-maxIdx") )
	)
	{
		usage (argv);
		return (0);
	}

	//filebase, min index, max index
	parse_argument (argc, argv, "-fb", fBase);
	parse_argument (argc, argv, "-minIdx", minIdx);
	parse_argument (argc, argv, "-maxIdx", maxIdx);

	//output transform file
	if (!find_switch(argc, argv, "-trFile") )
	{
		transformFileName = "trans.mat";
	}
	else
	{
		parse_argument (argc, argv, "-trFile", transformFileName);
	}

	//Skip count switch
	if (find_switch (argc, argv, "-s") )
	{
		parse_argument (argc, argv, "-s", stepSz);
		step = atoi(stepSz.c_str());
	}

	//transtype switch
	if (find_switch (argc, argv, "-t") )
	{
		hasOutfile = true;
		parse_argument (argc, argv, "-t", transStr);
		transType = getTranstype(transStr);
	}

	//gps kml switch, required for transType>=11
	if (find_switch (argc, argv, "-gps") )
	{
		hasGps = true;
		parse_argument (argc, argv, "-gps", gpsKmlFile);
		cout<<"Loading "<<gpsKmlFile<<endl;
		pf.loadKml(gpsKmlFile,"my path");
	}
	else if (transType >= 11)
	{
		cout << "Cannot run type GPS without kml gps input file <-gps [path-to-kml]>"<< endl;
		usage(argv);
		return(0);
	}


	//Intensity Filter switch
	if (find_switch (argc, argv, "-iFilt") )
	{
		parse_argument (argc, argv, "-iFilt", intStr);
		minInt = atoi(intStr.c_str());
		iFilter = true;
	}

	//load transform switch
	if (find_switch (argc, argv, "-lft") )
	{
		parse_argument (argc, argv, "-lft", ftFile);
		tmpTrans = tmpTrans*LoadFinalTransform(ftFile.c_str());
		cout << "Loaded starting transform:" << endl << tmpTrans << endl;
	}

	min = atoi(minIdx.c_str());
	max = atoi(maxIdx.c_str());

	START_TIMER(0);


	for (i=min;i<=max;i+=step)
	{
		START_TIMER(1);
		//cloud_ is the previous timestep
		//cloud2_ is the current timestep
		//swap
		cloud_->swap(*cloud2_);
		cloud2_->clear();
		//then load
		sprintf(fName,"%s/%04d.pcd",fBase.c_str(),i);
		pcl::io::loadPCDFile(fName,*cloud2_);

		//load the timestamp file, generated by PcapDump
		sprintf(fName,"%s/%04d.ts",fBase.c_str(),i);
		if ((tStampFile = fopen(fName,"r")))
		{
			fscanf(tStampFile,"%lld",&tStamp);
			if (i==min)startTime = tStamp / 1.0e6;

			deltaTime = (tStamp / 1.0e6) - (currTime + startTime);
			currTime = (tStamp / 1.0e6) - startTime;
			//cout << currTime << " " << deltaTime << " " << tStamp << endl;
			fclose(tStampFile);
		}
		else if (transType >= 12)
		{
			cout << fName << " not found, need timestamp for C-GPS\n" << endl;
			return 1;
		}

		sprintf(fName,"%s/%04d.pcd",fBase.c_str(),i);

    cout << "Processing " << fName << "..." << endl;

		if (iFilter)
		{
			filterIntensity(minInt,cloud2_);
		}

		if (!cloud_->empty())
		{
			tmpTransB = tmpTrans;

			//switch case for the tranformation types
			switch(transType)
			{
				case 0:
				//no tranform
				tmpTrans = Eigen::Matrix4f::Identity();
				break;
				case 1:
				//ICP trans
				tmpTrans = transform.getICPTransform(cloud_,cloud2_);
				break;
				case 2:
				//NDT trans
				tmpTrans = transform.getNDTTransform(cloud_,cloud2_);
				break;
				case 3:
				//PAIR trans
				tmpTrans = transform.getPairTransform(cloud_,cloud2_);
				break;
				case 4:
				//ADAPT trans
				tmpTrans = transform.getICPAdaptTransform(cloud_,cloud2_);
				break;
				case 5:
				//DoN trans
				tmpTrans = transform.getICPDoNTransform(cloud_,cloud2_);
				//tmpTrans = transform.getICPTransformGuess(cloud_,cloud2_,tmpTrans);
				break;
				case 6:
				//NDTDoN trans
				tmpTrans = transform.getNDTDoNTransform(cloud_,cloud2_);
				break;
				case 7:
				//Intensity filtered two pass ICP trans
				tmpTrans = transform.getICPIntTransform(cloud_,cloud2_);
				break;
				case 8:
				//Intensity filtered two pass ICP trans
				tmpTrans = transform.getICPVoxTransform(cloud_,cloud2_);
				break;
				case 9:
				tmpTrans = transform.getCubeTransform(cloud_,cloud2_,cSize,&cDist);
				break;
				case 10:
				//tmpTrans = transform.getCubeTransform(cloud_,cloud2_,tmpTrans,cDist,&cDist);
				//tmpTrans = transform.getICPTransformGuess(cloud_,cloud2_,tmpTrans,.05);
				tmpTrans = transform.getCubeICPTransform(cloud_,cloud2_,cSize,&cDist);
				break;
				case 11:
				tmpTrans = transform.getCubeGpsIcpTransform(cloud_,cloud2_,transMat,&pf,&cSize,&cDist,currTime-deltaTime,deltaTime);
				break;
				case 12:
				tmpTrans = transform.getCubeGpsTransform(cloud_,cloud2_,transMat,&pf,&cSize,&cDist,currTime-deltaTime,deltaTime);
				break;
				case 13:
				tmpTrans = transform.getGpsIcpTransform(cloud_,cloud2_,transMat,&pf,&cSize,currTime-deltaTime,deltaTime);
				break;
				case 14:
				tmpTrans = transform.getGpsIcpCubeTransform(cloud_,cloud2_,transMat,&pf,&cSize,&cDist,currTime-deltaTime,deltaTime);
				break;
				case 15:
				tmpTrans = transform.getGpsIcpNoCubeTransform(cloud_,cloud2_,transMat,&pf,&cSize,currTime-deltaTime,deltaTime);
				break;
				case 16:
				tmpTrans = transform.getGpsTransform(transMat,&pf,currTime-deltaTime,deltaTime);
				break;
				case 17:
				tmpTrans = transform.getCubeGpsIcpTransformAdapt(cloud_,cloud2_,transMat,&pf,currVel,currTime-deltaTime,deltaTime);
				break;
				default:
				//shouldn't happen
				break;
			}

		}

		//calculate distance traveled and acceleration
		prevVert = currVert;
		currVert = CVertex(transMat*tmpTrans);
		if (i > min+step)
		{
			prevVel = currVel;
			currVel = (currVert - prevVert);//*(1.0/deltaTime);
			if (i > min+(step*2))
			{
				currAcc = (currVel - prevVel);//*(1.0/deltaTime);
				//this could be attempted to watch for spikes in acceleration
				// to determine a poor alignment and try a different alignment technique
				/*if (currAcc.getDist() > 10.0)// && transType == getTranstype(transStr))
				{
				//transType = 11;
				cloud2_->swap(*cloud_);
				currTime -= deltaTime;
				currVert = prevVert;
				currVel = prevVel;
				//cout << "0 0 0 0 0 0 0 0"<< endl;
				cerr << "Skipping file " << i << endl;
				continue;
			}*/
		}
	}
	// if transtype is switched above, it needs to be reset
	//transType = getTranstype(transStr);

	//print out the distance and gps expected distance for the timestep
	// data available above to print acceleration magnitued with currAcc.getDist();
	if ( hasGps )
	cout << i << " " <<  currVel.getDist() << " " << pf.getDeltaForTime((currTime-deltaTime),deltaTime).norm() << endl;
	else
	cout << i << " " << currVel.getDist() << " processing...";


	//Use gps to rotate direction to match path (poor results)
	if (hasGps && transType < 11)
	{
		if (i>min)
		{
			prevPoint = pf.extractOrigin(transMat*tmpTransB);
			currPoint = pf.extractOrigin(transMat*tmpTransB*tmpTrans);
			if (!pf.isInSegment(segIdx,currPoint))
			{
				segIdx++;
				cout << "Next Segment after " << (i-min) << " points" << endl;
			}
			tmpTransB = tmpTransB*pf.rotateAlignment(prevPoint,currPoint,segIdx);
			transMat = transMat*tmpTransB;
			sprintf(fNameB,"%s/%04d.pcd",fBase.c_str(),i-step);
			SaveTransform(fNameB,tmpTransB);
		}
	}
	//the first transform with a gps file alignment has to be rotated
	// so the tranformations start following the path of the gps
	// this allows the estimated velocities to be in the correct directions
	else if (hasGps && i==min)
	{
		PointXYZ p1,p2;
		p1.x=0;p1.y=0;p1.z=0;
		p2.x=0;p2.y=1;p2.z=0;
		firstRot = pf.rotateAlignment2D(p1,p2,0);
		transMat = transMat*firstRot;
		sprintf(fNameB,"%s/%04d.pcd",fBase.c_str(),i);
		SaveTransform(fNameB,firstRot);
	}
	else
	{
		transMat = transMat*tmpTrans;
		SaveTransform(fName,tmpTrans);
	}

	//Timer in fractional seconds for length of 1 alignment
	STOP_TIMER(1,timerString);
}
//Timer in fractional seconds for length of all alignments
STOP_TIMER(0,longTimerString);
cout << longTimerString << "\n";

//in this case we still need to save the last result
if (hasGps && transType < 11)
{
	SaveTransform(fName,tmpTrans);
}

//save out the final transform for resuming past the origin
// loading won't work with a gps path, it could if the start and curr times
// were saved and used to determine the gps velocity vector when starting
// again
if (find_switch (argc, argv, "-sft") )
{
	parse_argument (argc, argv, "-sft", ftFile);
	SaveFinalTransform(ftFile.c_str(),transMat*tmpTrans);
}

fclose(transFile);

return (0);
}
