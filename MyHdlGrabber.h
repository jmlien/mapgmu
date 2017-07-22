/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2012 The MITRE Corporation
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
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVERsignal
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "pcl/pcl_config.h"

#ifndef MY_HDL_GRABBER_H_
#define MY_HDL_GRABBER_H_

#include <pcl/io/grabber.h>
#include <pcl/io/impl/synchronized_queue.hpp>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <boost/asio.hpp>
#include <string>
#include <pcap/pcap.h>


  class MyHdlGrabber : public pcl::Grabber
  {
    public:
      
      MyHdlGrabber (const std::string& correctionsFile = "",
	const std::string& pcapFile = "");

      /** \brief virtual Destructor inherited from the Grabber interface. It never throws. */
      virtual ~MyHdlGrabber () throw ();
      
      /** \brief Starts processing the Velodyne packets, either from the network or PCAP file. */
      virtual void start ();

      /** \brief Stops processing the Velodyne packets, either from the network or PCAP file */
      virtual void stop ();
      virtual bool isRunning() const;
      virtual std::string getName() const;
      virtual float getFramesPerSecond() const;
      bool hasClouds();
      pcl::PointCloud<pcl::PointXYZI>::Ptr getNextCloud();

    protected:
      static const int HDL_NUM_ROT_ANGLES = 36001;
      static const int HDL_LASER_PER_FIRING = 32;
      static const int HDL_MAX_NUM_LASERS = 64;
      static const int HDL_FIRING_PER_PKT = 12;
      
      
      enum HDLBlock
       {
         BLOCK_0_TO_31 = 0xeeff, BLOCK_32_TO_63 = 0xddff
       };
      
 #pragma pack(push, 1)
       typedef struct HDLLaserReturn
       {
           unsigned short distance;
           unsigned char intensity;
       } HDLLaserReturn;
 #pragma pack(pop)
 
       struct HDLFiringData
       {
           unsigned short blockIdentifier;
           unsigned short rotationalPosition;
           HDLLaserReturn laserReturns[HDL_LASER_PER_FIRING];
       };
 
       struct HDLDataPacket
       {
           HDLFiringData firingData[HDL_FIRING_PER_PKT];
           unsigned int gpsTimestamp;
           unsigned char blank1;
           unsigned char blank2;
       };
 
       struct HDLLaserCorrection
       {
           double azimuthCorrection;
           double verticalCorrection;
           double distanceCorrection;
           double verticalOffsetCorrection;
           double horizontalOffsetCorrection;
           double sinVertCorrection;
           double cosVertCorrection;
           double sinVertOffsetCorrection;
           double cosVertOffsetCorrection;
       };

    private:
      struct pcap_pkthdr *header;
      const unsigned char *data;
      pcap_t *pcapPtr;
      int sweepCount;
      int cloudCount;
      int lastOffset;
      bool with_image;
      static double *cos_lookup_table_;
      static double *sin_lookup_table_;
      boost::asio::ip::udp::endpoint udp_listener_endpoint_;
      boost::asio::ip::address source_address_filter_;
      unsigned short source_port_filter_;
	    boost::asio::io_service hdl_read_socket_service_;
      boost::asio::ip::udp::socket *hdl_read_socket_;
      std::string my_pcap_file_;
      boost::thread *pcap_thread;
      HDLLaserCorrection laser_corrections_[HDL_MAX_NUM_LASERS];
      bool terminate_;
      bool running_;
      int last_azimuth_;
      float min_distance_threshold_;
      float max_distance_threshold_;
      
      
      void initialize (const std::string& correctionsFile);
      void loadCorrectionsFile (const std::string& correctionsFile);
      void loadHDL32Corrections ();

      void pcapThread();
      bool completePointCloud (HDLDataPacket *dataPacket,
			       pcl::PointCloud<pcl::PointXYZI>::Ptr myCloud);
      void computeXYZI (pcl::PointXYZI& pointXYZI, int azimuth,
          HDLLaserReturn laserReturn, HDLLaserCorrection correction);
  };

#endif 
