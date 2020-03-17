//
// Created by tlb on 18-9-11.
//
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl/io/pcd_io.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <tf/tf.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

#include "backward.hpp"

using namespace std;
using namespace Eigen;
ros::Subscriber _map_sub;

void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);
int i = 0;
int lastsize = 0;
int main(int argc, char **argv) {
  ros::init(argc, argv, "map_recorder");
  ros::NodeHandle nh("~");
  _map_sub = nh.subscribe("/random_map_sensing/all_map", 1, rcvPointCloudCallBack);
  while (ros::ok())
    ros::spinOnce();
  return 0;
}
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map) {
  pcl::PointCloud<pcl::PointXYZ> cloud;
  pcl::fromROSMsg(pointcloud_map, cloud);

  if ((int) cloud.points.size() == 0)
    return;
  stringstream dir;
  dir << "/home/tlb/catkin_test/src/Btraj/data/map" << i + 1 << ".pcd";
  if (cloud.points.size() == lastsize)
    return;
  lastsize = cloud.points.size();
  pcl::io::savePCDFileASCII(dir.str(), cloud);
  std::cerr << "Saved the " << ++i << "maps" << cloud.points.size() << " data points to map.pcd." << std::endl;

}