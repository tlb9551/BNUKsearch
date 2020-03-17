#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>

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

#include "bnuk_search.h"
#include "b_spline_optimizer.h"
#include "uniform_bspline_3d.h"
#include "quadrotor_msgs/BSplineTrajectory.h"
#include "quadrotor_msgs/PositionCommand.h"

#include <chrono>

using namespace std;
using namespace Eigen;
using namespace sdf_tools;
using std::chrono::high_resolution_clock;

namespace backward {
backward::SignalHandling sh;
}

// simulation param from launch file
bool _use_opt;
double _vis_traj_width;
double _resolution, _inv_resolution;
double _cloud_margin, _cube_margin, _check_horizon, _stop_horizon;//_cloud_margin 点云膨胀系数（点云膨胀到周围多少）
double _x_size, _y_size, _z_size, _x_local_size, _y_local_size, _z_local_size;
double _localmap_x_size, _localmap_y_size, _localmap_z_size;
Vector3d end_pt_global;
Vector3d end_vel_global;
Vector3d _current_traj_start_pt;

///bnuk_log
struct bnuk_log_datatype {
  vector<Eigen::Vector3d> traj_pos;
  vector<double> traj_time;
  vector<double> traj_lenth_vec;
  double traj_lenth;
  double traj_duration;

  vector<double> compute_time;
  double compute_time_mean;
  int replan;
  void clear() {
    traj_time.clear();
    traj_pos.clear();
    compute_time.clear();
    traj_lenth_vec.clear();

    traj_lenth = 0.0;
    traj_duration = 0.0;
    compute_time_mean = 0.0;
  }
} bnuk_log, bnuk_opt_log;
//////BNUK
vector<Vector3d> _start_state;
vector<Vector3d> _end_state;
double _tie_breaker;
double _heu_coeff;
double _max_vel(1.6);
double _dynamics_max_vel(1.6), _dynamics_max_acc(1.6);
double cw_jerk(0.004);
double _delt_t(1.0);
ewok::UniformBSpline3D<6, double>::Ptr _b_spline_, _b_spline_opt_;
shared_ptr<bnuk_search> _searcher;
bool kinodynamic_traj_planning();
void visBSplineTrajectory(ewok::UniformBSpline3D<6, double>::Ptr &b_spline_, string ns);
void visBSplineControlPts(ewok::UniformBSpline3D<6, double>::Ptr &b_spline_, string ns);
void visExpNode(vector<BSplineGridNodePtr> nodes);
void visVel();
void bnuk_log_write_VFRBK(const std::string file_name);
bool solve_opt();

bool _is_use_fm, _is_proj_cube, _is_limit_vel, _is_limit_acc;
int _max_inflate_iter, _traj_order;
double _minimize_order;

// useful global variables
nav_msgs::Odometry _odom;
bool _has_odom = false;
bool _has_map = false;
bool _has_target = false;
bool _has_traj = false;
bool _is_emerg = false;
bool _is_init = true;

Vector3d _start_pt, _start_vel, _start_acc, _end_pt, _end_vel;

double _init_x, _init_y, _init_z;
Vector3d _map_origin;
double _pt_max_x, _pt_min_x, _pt_max_y, _pt_min_y, _pt_max_z, _pt_min_z;
int _max_x_id, _max_y_id, _max_z_id, _max_local_x_id, _max_local_y_id, _max_local_z_id;
int _traj_id = 1;
COLLISION_CELL _free_cell(0.0);
COLLISION_CELL _obst_cell(1.0);
// ros related
ros::Subscriber _map_sub, _pts_sub, _odom_sub;
ros::Publisher _local_map_vis_pub, _inf_map_vis_pub, _corridor_vis_pub, _traj_vis_pub, _grid_path_vis_pub,
    _nodes_vis_pub, _traj_pub, _checkTraj_vis_pub, _stopTraj_vis_pub, _end_vis_pub;

quadrotor_msgs::BSplineTrajectory _traj;
ros::Time _start_time = ros::TIME_MAX;
shared_ptr<CollisionMapGrid> collision_map_local;

void rcvWaypointsCallback(const nav_msgs::Path &wp);
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map);
void rcvOdometryCallbck(const nav_msgs::Odometry odom);
void setMidEndState();

bool checkExecTraj();
bool checkCoordObs(Vector3d checkPt);

quadrotor_msgs::BSplineTrajectory getBsplineTraj(ewok::UniformBSpline3D<6, double>::Ptr bspline);

int odom_count(0);
ros::Time last_odom;
Vector3d traj_pre;
void rcvOdometryCallbck(const nav_msgs::Odometry odom) {
  if (odom.pose.pose.position.x != -20.0 && odom.pose.pose.position.y != -20.0) {
    if (odom_count >= 4) {

      if (last_odom.toSec() != 0.0) {
        //cout<<"record"<<endl;
        odom_count = 0;
        Vector3d traj_cur;
        traj_cur(0) = _odom.pose.pose.position.x;
        traj_cur(1) = _odom.pose.pose.position.y;
        traj_cur(2) = _odom.pose.pose.position.z;
        bnuk_log.traj_pos.push_back(traj_cur);

        double length_inc = (traj_cur - traj_pre).norm();
        bnuk_log.traj_lenth += length_inc;
        bnuk_log.traj_lenth_vec.push_back(bnuk_log.traj_lenth);

        if (length_inc > 0.0001) {
          bnuk_log.traj_duration += (_odom.header.stamp - last_odom).toSec();
          bnuk_log.traj_time.push_back(bnuk_log.traj_duration);
        }

        last_odom = _odom.header.stamp;
        traj_pre = traj_cur;

      } else {
        last_odom = _odom.header.stamp;
        traj_pre(0) = _odom.pose.pose.position.x;
        traj_pre(1) = _odom.pose.pose.position.y;
        traj_pre(2) = _odom.pose.pose.position.z;
      }
    } else
      odom_count++;
  }

  if (odom.header.frame_id != "uav")
    return;

  _odom = odom;
  _has_odom = true;

  _start_pt(0) = _odom.pose.pose.position.x;
  _start_pt(1) = _odom.pose.pose.position.y;
  _start_pt(2) = _odom.pose.pose.position.z;

  _start_vel(0) = _odom.twist.twist.linear.x;
  _start_vel(1) = _odom.twist.twist.linear.y;
  _start_vel(2) = _odom.twist.twist.linear.z;

  _start_acc(0) = _odom.twist.twist.angular.x;
  _start_acc(1) = _odom.twist.twist.angular.y;
  _start_acc(2) = _odom.twist.twist.angular.z;

  if (std::isnan(_odom.pose.pose.position.x) || std::isnan(_odom.pose.pose.position.y)
      || std::isnan(_odom.pose.pose.position.z))
    return;

  static tf::TransformBroadcaster br;
  tf::Transform transform;
  transform.setOrigin(tf::Vector3(_odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z));
  transform.setRotation(tf::Quaternion(0, 0, 0, 1.0));
  br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "quadrotor"));
}

void rcvWaypointsCallback(const nav_msgs::Path &wp) {
  if (wp.poses[0].pose.position.z < 0.0)
    return;

  end_pt_global << wp.poses[0].pose.position.x,
      wp.poses[0].pose.position.y,
      wp.poses[0].pose.position.z;
//    end_vel_global <<2.0*wp.poses[0].pose.orientation.w*wp.poses[0].pose.orientation.w-1,
//                    2.0*wp.poses[0].pose.orientation.z*wp.poses[0].pose.orientation.w,
//                    0.0;
  end_vel_global << 0.0, 0.0, 0.0;
  Vector3d local_map_origin = collision_map_local->GetOriginTransform().translation();
  cout << "localmap_origin" << local_map_origin.transpose() << endl;
  Vector3d start = collision_map_local->GetOriginTransform().translation()
      + Vector3d(_localmap_x_size / 2.0, _localmap_y_size / 2.0, _localmap_z_size / 2.0);
  if (end_pt_global(0) >= local_map_origin(0) && end_pt_global(0) <= local_map_origin(0) + _localmap_x_size
      && end_pt_global(1) >= local_map_origin(1) && end_pt_global(1) <= local_map_origin(1) + _localmap_y_size
      && end_pt_global(2) >= local_map_origin(2)
      && end_pt_global(2) <= local_map_origin(2) + _localmap_z_size) {//target in local map
    cout << "target in local map " << end_pt_global.transpose() << endl;
    _end_pt = end_pt_global;
    _end_vel = Vector3d::Zero();
  } else {    //target not in local map
    cout << "target not in local map " << end_pt_global.transpose() << endl;

    Vector3d d;
    d = end_pt_global - start;
    double k1(0.0), k2(0.0), k3(0.0), k4(0.0), k5(0.0), k6(0.0);
    if (d(0) != 0.0) {
      k1 = (local_map_origin(0) + _localmap_x_size - start(0)) / d(0);
      k4 = (local_map_origin(0) - start(0)) / d(0);
    } else {
      k1 = 10;
      k4 = 10;
    }
    if (d(1) != 0.0) {
      k2 = (local_map_origin(1) + _localmap_y_size - start(1)) / d(1);
      k5 = (local_map_origin(1) - start(1)) / d(1);
    } else {
      k2 = 10;
      k5 = 10;
    }
    if (d(2) != 0.0) {
      k3 = (local_map_origin(2) + _localmap_z_size - start(2)) / d(2);
      k6 = (local_map_origin(2) - start(2)) / d(2);
    } else {
      k3 = 10;
      k6 = 10;
    }

    k1 = k1 > 0 ? k1 : 10;
    k2 = k2 > 0 ? k2 : 10;
    k3 = k3 > 0 ? k3 : 10;
    k4 = k4 > 0 ? k4 : 10;
    k5 = k5 > 0 ? k5 : 10;
    k6 = k6 > 0 ? k6 : 10;

    double k = min(k1, min(k2, min(k3, min(k4, min(k5, k6)))));
    Vector3d insect_point = start + k * (end_pt_global - start) * 0.95;
    Vector3d rest_d = (1.0 - k) * (end_pt_global - start);
    Vector3d velocity;

    _end_pt = insect_point;
    velocity = rest_d.normalized() * _max_vel / 3.0;
    _end_vel = velocity;
  }

  _has_target = true;
  _is_emerg = true;

  ROS_INFO("[BNUK Node] receive the way-points");

  while (!kinodynamic_traj_planning()) {
    ros::Duration(0.05).sleep();
  };
}
void setMidEndState() {

  Vector3d local_map_origin = collision_map_local->GetOriginTransform().translation();
  cout << "localmap_origin" << local_map_origin.transpose() << endl;
  Vector3d start = collision_map_local->GetOriginTransform().translation()
      + Vector3d(_localmap_x_size / 2.0, _localmap_y_size / 2.0, _localmap_z_size / 2.0);
  if (end_pt_global(0) >= local_map_origin(0) && end_pt_global(0) <= local_map_origin(0) + _localmap_x_size
      && end_pt_global(1) >= local_map_origin(1) && end_pt_global(1) <= local_map_origin(1) + _localmap_y_size
      && end_pt_global(2) >= local_map_origin(2)
      && end_pt_global(2) <= local_map_origin(2) + _localmap_z_size) {//target in local map
    cout << "target in local map " << end_pt_global.transpose() << endl;
    _end_pt = end_pt_global;
    _end_vel = Vector3d::Zero();
  } else {    //target not in local map
    cout << "target not in local map " << end_pt_global.transpose() << endl;

    Vector3d d;
    d = end_pt_global - start;
    double k1(0.0), k2(0.0), k3(0.0), k4(0.0), k5(0.0), k6(0.0);
    if (d(0) != 0.0) {
      k1 = (local_map_origin(0) + _localmap_x_size - start(0)) / d(0);
      k4 = (local_map_origin(0) - start(0)) / d(0);
    } else {
      k1 = 10;
      k4 = 10;
    }
    if (d(1) != 0.0) {
      k2 = (local_map_origin(1) + _localmap_y_size - start(1)) / d(1);
      k5 = (local_map_origin(1) - start(1)) / d(1);
    } else {
      k2 = 10;
      k5 = 10;
    }
    if (d(2) != 0.0) {
      k3 = (local_map_origin(2) + _localmap_z_size - start(2)) / d(2);
      k6 = (local_map_origin(2) - start(2)) / d(2);
    } else {
      k3 = 10;
      k6 = 10;
    }

    k1 = k1 > 0 ? k1 : 10;
    k2 = k2 > 0 ? k2 : 10;
    k3 = k3 > 0 ? k3 : 10;
    k4 = k4 > 0 ? k4 : 10;
    k5 = k5 > 0 ? k5 : 10;
    k6 = k6 > 0 ? k6 : 10;

    double k = min(k1, min(k2, min(k3, min(k4, min(k5, k6)))));
    Vector3d insect_point = start + k * (end_pt_global - start) * 0.95;
    Vector3d rest_d = (1.0 - k) * (end_pt_global - start);
    Vector3d velocity;

    _end_pt = insect_point;
    velocity = rest_d.normalized() * _max_vel / 3.0;
    _end_vel = velocity;
  }

  _has_target = true;
  _is_emerg = true;

  ROS_INFO("[BNUK Node] mid  way-points");
}

Vector3d _local_origin;
/*接受点云，膨胀点云，建立障碍物地图，trajPlanning（）*/
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 &pointcloud_map) {
  pcl::PointCloud<pcl::PointXYZ> cloud;
  pcl::fromROSMsg(pointcloud_map, cloud);

  if ((int) cloud.points.size() == 0)
    return;

  double local_c_x = (int) ((_start_pt(0) - _x_local_size / 2.0) * _inv_resolution + 0.5) * _resolution;
  double local_c_y = (int) ((_start_pt(1) - _y_local_size / 2.0) * _inv_resolution + 0.5) * _resolution;
  double local_c_z = (int) ((_start_pt(2) - _z_local_size / 2.0) * _inv_resolution + 0.5) * _resolution;

  _local_origin << local_c_x, local_c_y, local_c_z;

  Translation3d origin_local_translation(_local_origin(0), _local_origin(1), _local_origin(2));
  Quaterniond origin_local_rotation(1.0, 0.0, 0.0, 0.0);

  Affine3d origin_local_transform = origin_local_translation * origin_local_rotation;

  collision_map_local = make_shared<CollisionMapGrid>(origin_local_transform,
                                                      "world",
                                                      _resolution,
                                                      _localmap_x_size,
                                                      _localmap_y_size,
                                                      _localmap_z_size,
                                                      _free_cell);

  vector<pcl::PointXYZ> inflatePts(20);
  pcl::PointCloud<pcl::PointXYZ> cloud_inflation;
  pcl::PointCloud<pcl::PointXYZ> cloud_local;

  for (int idx = 0; idx < (int) cloud.points.size(); idx++) {
    auto mk = cloud.points[idx];
    pcl::PointXYZ pt(mk.x, mk.y, mk.z);

    if (fabs(pt.x - _start_pt(0)) > _x_local_size / 2.0 || fabs(pt.y - _start_pt(1)) > _y_local_size / 2.0
        || fabs(pt.z - _start_pt(2)) > _z_local_size / 2.0)
      continue;

    cloud_local.push_back(pt);
    Vector3d addPt(pt.x, pt.y, pt.z);
    collision_map_local->Set3d(addPt, _obst_cell);
  }

  if (_start_pt(0) - _localmap_x_size / 2.0 < -_x_size / 2.0) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(-_x_size / 2.0 + 0.5 * _resolution, 0.0, 0.0));
    for (int64_t y = 0; y < collision_map_local->GetNumYCells(); y++)
      for (int64_t z = 0; z < collision_map_local->GetNumZCells(); z++)
        collision_map_local->Set(int64_t(bound(0)), y, z, _obst_cell);
  }
  if (_start_pt(0) + _localmap_x_size / 2.0 > _x_size / 2.0) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(_x_size / 2.0 - 0.5 * _resolution, 0.0, 0.0));
    for (int64_t y = 0; y < collision_map_local->GetNumYCells(); y++)
      for (int64_t z = 0; z < collision_map_local->GetNumZCells(); z++)
        collision_map_local->Set(int64_t(bound(0)), y, z, _obst_cell);
  }

  if (_start_pt(1) - _localmap_y_size / 2.0 < -_y_size / 2.0) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(0.0, -_y_size / 2.0 + 0.5 * _resolution, 0.0));
    for (int64_t x = 0; x < collision_map_local->GetNumXCells(); x++)
      for (int64_t z = 0; z < collision_map_local->GetNumZCells(); z++)
        collision_map_local->Set(x, int64_t(bound(1)), z, _obst_cell);
  }
  if (_start_pt(1) + _localmap_y_size / 2.0 > _y_size / 2.0) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(0.0, _y_size / 2.0 - 0.5 * _resolution, 0.0));
    for (int64_t x = 0; x < collision_map_local->GetNumXCells(); x++)
      for (int64_t z = 0; z < collision_map_local->GetNumZCells(); z++)
        collision_map_local->Set(x, int64_t(bound(1)), z, _obst_cell);
  }

  if (_start_pt(2) - _localmap_z_size / 2.0 < 0.0) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(0.0, 0.0, 0.5 * _resolution));
    for (int64_t x = 0; x < collision_map_local->GetNumXCells(); x++)
      for (int64_t y = 0; y < collision_map_local->GetNumYCells(); y++)
        collision_map_local->Set(x, y, int64_t(bound(2)), _obst_cell);
  }
  if (_start_pt(2) + _localmap_z_size / 2.0 > _z_size) {
    Vector3i bound = collision_map_local->LocationToGridIndex(Vector3d(0.0, 0.0, _z_size - 0.5 * _resolution));
    for (int64_t x = 0; x < collision_map_local->GetNumXCells(); x++)
      for (int64_t y = 0; y < collision_map_local->GetNumYCells(); y++)
        collision_map_local->Set(x, y, int64_t(bound(2)), _obst_cell);
  }

  _has_map = true;

  cloud_local.width = cloud_local.points.size();
  cloud_local.height = 1;
  cloud_local.is_dense = true;
  cloud_local.header.frame_id = "world";

  sensor_msgs::PointCloud2 localMap;
  pcl::toROSMsg(cloud_local, localMap);
  _local_map_vis_pub.publish(localMap);

  if (checkExecTraj() == true)//前面有碰撞
    while (!kinodynamic_traj_planning()) {
      ros::Duration(0.05).sleep();
    };
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "bnuk_node");
  ros::NodeHandle nh("~");

  _map_sub = nh.subscribe("/random_forest_sensing/random_forest", 1, rcvPointCloudCallBack);
  _odom_sub = nh.subscribe("/odom/fake_odom", 1, rcvOdometryCallbck);
  _pts_sub = nh.subscribe("/waypoint_generator/waypoints", 1, rcvWaypointsCallback);

  _inf_map_vis_pub = nh.advertise<sensor_msgs::PointCloud2>("vis_map_inflate", 1);
  _local_map_vis_pub = nh.advertise<sensor_msgs::PointCloud2>("vis_map_local", 1);
  _traj_vis_pub = nh.advertise<visualization_msgs::Marker>("trajectory_vis", 1);
  _corridor_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("corridor_vis", 1);
  _grid_path_vis_pub = nh.advertise<visualization_msgs::MarkerArray>("grid_path_vis", 1);
  _nodes_vis_pub = nh.advertise<visualization_msgs::Marker>("expanded_nodes_vis", 1);
  _checkTraj_vis_pub = nh.advertise<visualization_msgs::Marker>("check_trajectory", 1);
  _stopTraj_vis_pub = nh.advertise<visualization_msgs::Marker>("stop_trajectory", 1);
  _end_vis_pub = nh.advertise<visualization_msgs::Marker>("end_vector", 1);

  _traj_pub = nh.advertise<quadrotor_msgs::BSplineTrajectory>("trajectory", 10);

  nh.param("map/resolution", _resolution, 0.2);

  nh.param("map/x_size", _x_size, 50.0);
  nh.param("map/y_size", _y_size, 50.0);
  nh.param("map/z_size", _z_size, 5.0);

  nh.param("map/x_local_size", _x_local_size, 20.0);
  nh.param("map/y_local_size", _y_local_size, 20.0);
  nh.param("map/z_local_size", _z_local_size, 5.0);

  nh.param("planning/init_x", _init_x, -20.0);
  nh.param("planning/init_y", _init_y, -20.0);
  nh.param("planning/init_z", _init_z, 0.0);
  nh.param("planning/tie_breaker", _tie_breaker, 1.03);

  nh.param("planning/max_vel", _max_vel, 1.6);
  nh.param("planning/dynamic_max_acc", _dynamics_max_vel, 1.6);
  nh.param("planning/dynamic_max_vel", _dynamics_max_acc, 1.6);

  nh.param("planning/check_horizon", _check_horizon, 10.0);
  nh.param("planning/stop_horizon", _stop_horizon, 3.0);

  nh.param("vis/vis_traj_width", _vis_traj_width, 0.15);
  nh.param("use_opt", _use_opt, false);

  _b_spline_.reset(new ewok::UniformBSpline3D<6, double>(_delt_t));
  _b_spline_opt_.reset(new ewok::UniformBSpline3D<6, double>(_delt_t));
  _map_origin << -_x_size / 2.0, -_y_size / 2.0, 0.0;
  _pt_max_x = +_x_size / 2.0;
  _pt_min_x = -_x_size / 2.0;
  _pt_max_y = +_y_size / 2.0;
  _pt_min_y = -_y_size / 2.0;
  _pt_max_z = +_z_size;
  _pt_min_z = 0.0;

  _inv_resolution = 1.0 / _resolution;
  _max_x_id = (int) (_x_size * _inv_resolution);
  _max_y_id = (int) (_y_size * _inv_resolution);
  _max_z_id = (int) (_z_size * _inv_resolution);
  _max_local_x_id = (int) (_x_local_size * _inv_resolution);
  _max_local_y_id = (int) (_y_local_size * _inv_resolution);
  _max_local_z_id = (int) (_z_local_size * _inv_resolution);

  double _buffer_size = 1.6;
  _localmap_x_size = _x_local_size + _buffer_size;
  _localmap_y_size = _y_local_size + _buffer_size;
  _localmap_z_size = _z_local_size;

  Translation3d origin_translation(_map_origin(0), _map_origin(1), 0.0);
  Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
  Affine3d origin_transform = origin_translation * origin_rotation;
  _searcher = make_shared<bnuk_search>(Vector3i(_localmap_x_size * _inv_resolution,
                                                _localmap_y_size * _inv_resolution,
                                                _localmap_z_size * _inv_resolution), 1.0);
  ros::Rate rate(100);
  bool status = ros::ok();
  while (status) {
    nh.getParam("planning/tie_breaker", _tie_breaker);
    bool log;
    nh.getParam("log", log);
    if (log)
      bnuk_log_write_VFRBK("/home/tlb/catkin_test/src/BNUK_Opt/comparison/replan/bnuk.txt");

    ros::spinOnce();
    status = ros::ok();
    rate.sleep();
  }

  return 0;
}

bool checkExecTraj() {
  if (_has_traj == false)
    return false;

  if (abs((_start_pt - _current_traj_start_pt).maxCoeff()) > _x_local_size / 2.0 - 2 * _max_vel) {
    setMidEndState();
    return true;
  }

  Vector3d traj_pt;

  visualization_msgs::Marker _check_traj_vis, _stop_traj_vis;

  geometry_msgs::Point pt;
  _stop_traj_vis.header.stamp = _check_traj_vis.header.stamp = ros::Time::now();
  _stop_traj_vis.header.frame_id = _check_traj_vis.header.frame_id = "world";

  _check_traj_vis.ns = "trajectory/check_trajectory";
  _stop_traj_vis.ns = "trajectory/stop_trajectory";

  _stop_traj_vis.id = _check_traj_vis.id = 0;
  _stop_traj_vis.type = _check_traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
  _stop_traj_vis.action = _check_traj_vis.action = visualization_msgs::Marker::ADD;

  _stop_traj_vis.scale.x = 2.0 * _vis_traj_width;
  _stop_traj_vis.scale.y = 2.0 * _vis_traj_width;
  _stop_traj_vis.scale.z = 2.0 * _vis_traj_width;

  _check_traj_vis.scale.x = 1.5 * _vis_traj_width;
  _check_traj_vis.scale.y = 1.5 * _vis_traj_width;
  _check_traj_vis.scale.z = 1.5 * _vis_traj_width;

  _check_traj_vis.pose.orientation.x = 0.0;
  _check_traj_vis.pose.orientation.y = 0.0;
  _check_traj_vis.pose.orientation.z = 0.0;
  _check_traj_vis.pose.orientation.w = 1.0;

  _stop_traj_vis.pose = _check_traj_vis.pose;

  _stop_traj_vis.color.r = 0.0;
  _stop_traj_vis.color.g = 1.0;
  _stop_traj_vis.color.b = 0.0;
  _stop_traj_vis.color.a = 1.0;

  _check_traj_vis.color.r = 0.0;
  _check_traj_vis.color.g = 0.0;
  _check_traj_vis.color.b = 1.0;
  _check_traj_vis.color.a = 1.0;
  double t_now = max(0.0, (_odom.header.stamp - _start_time).toSec());
  double max_valid_time;
  if (_use_opt) {
    max_valid_time = _b_spline_opt_->maxValidTime();
  } else {
    max_valid_time = _b_spline_->maxValidTime();
  }
  for (double t = t_now; t < _b_spline_->maxValidTime(); t += 0.01) {
    if (t > _check_horizon + t_now) break;
    if (_use_opt) {
      traj_pt = _b_spline_opt_->evaluate(t, 0);

    } else {
      traj_pt = _b_spline_->evaluate(t, 0);
    }
    pt.x = traj_pt(0);
    pt.y = traj_pt(1);
    pt.z = traj_pt(2);
    //pt odom时间bezier轨迹曲线上点的pose
    _check_traj_vis.points.push_back(pt);

    if (t <= _stop_horizon + t_now)
      _stop_traj_vis.points.push_back(pt);

    if (checkCoordObs(traj_pt))//检测到障碍物在traj_pt点
    {
      ROS_WARN("predicted collision time is %f ahead", t);

      if (t <= _stop_horizon + t_now) {
        ROS_ERROR("emergency occurs in time is %f ahead", t);
        _is_emerg = true;
      }

      _checkTraj_vis_pub.publish(_check_traj_vis);
      _stopTraj_vis_pub.publish(_stop_traj_vis);

      setMidEndState();
      return true;
    }
  }
  _checkTraj_vis_pub.publish(_check_traj_vis);
  _stopTraj_vis_pub.publish(_stop_traj_vis);

  return false;
}
bool checkCoordObs(Vector3d checkPt) {
  if (collision_map_local->Get(checkPt(0), checkPt(1), checkPt(2)).first.occupancy > 0.5)
    return true;

  return false;
}

bool kinodynamic_traj_planning() {
  ros::Time time_vf1 = ros::Time::now();

  _searcher->initGridNodeMap(_resolution);

  if (_has_target == false || _has_map == false || _has_odom == false)
    return false;


  /*** collision map A* ->  control points**/
  /**control point -> b-spline trajectory***/
  _b_spline_->ClearControlPoints();
  _b_spline_opt_->ClearControlPoints();
  _searcher->setTieBreaker(_tie_breaker);
  _searcher->setHeuCoeff(_heu_coeff);
  _searcher->setMaxVel(_max_vel);
  _searcher->setDynamicsConstraint(Vector4f(_dynamics_max_vel,
                                            -_dynamics_max_vel,
                                            _dynamics_max_acc,
                                            -_dynamics_max_acc));
  _searcher->setQuadraticCostWeight(std::array<double, 4>{0.0, 0.0, cw_jerk, 0.0});
  _searcher->linkCollsionMap(collision_map_local);
  ros::Time time_vf2 = ros::Time::now();
  _start_state.clear();
  _end_state.clear();
  _start_state.push_back(_start_pt);
  _start_state.push_back(_start_vel);
  _end_state.push_back(_end_pt);
  _end_state.push_back(_end_vel);
  _searcher->setStartState(_start_state);
  _searcher->setEndState(_end_state);
  //visVel();
  _searcher->V_F_RBK_search();
  ros::Time time_search2 = ros::Time::now();
  bnuk_log.compute_time.push_back((time_search2 - time_vf1).toNSec() / 1000000.0);
  cout << "bnuk_search compute_time:" << (time_search2 - time_vf1).toNSec() / 1000000.0 << "ms" << endl;
  bnuk_log.compute_time_mean =
      accumulate(bnuk_log.compute_time.begin(), bnuk_log.compute_time.end(), 0.0) / bnuk_log.compute_time.size();
  bnuk_log.replan = bnuk_log.compute_time.size() - 1;
  vector<Eigen::Vector3d> ControlPts = _searcher->getpath();
  cout << "controlpoints num " << ControlPts.size() << endl;
  for (auto it = ControlPts.begin(); it != ControlPts.end(); ++it) {
    _b_spline_->push_back(*it);
  }
  bool res_opt = false;
  if (_use_opt) {
    res_opt = solve_opt();
    if (res_opt) {
      _is_emerg = false;
      _has_traj = true;
    } else {
      _is_emerg = false;
      _has_traj = false;
      return false;
    }
  } else {
    _is_emerg = false;
    _has_traj = true;
  }

  if (_use_opt && res_opt) {
    visBSplineControlPts(_b_spline_opt_, "bnuk_opt");
    visBSplineTrajectory(_b_spline_opt_, "bnuk_opt");
    _traj = getBsplineTraj(_b_spline_opt_);
    _traj_pub.publish(_traj);
    _current_traj_start_pt = _start_pt;
    _traj_id++;
  } else {
    visBSplineControlPts(_b_spline_, "bnuk");
    visBSplineTrajectory(_b_spline_, "bnuk");
    _traj = getBsplineTraj(_b_spline_);
    _traj_pub.publish(_traj);
    _current_traj_start_pt = _start_pt;
    _traj_id++;
  }

  /****visExpNode***/
  auto expnodes = _searcher->getVisitedNodes();
  visExpNode(expnodes);


  /****visVelocityField***/

  _searcher->resetLocalMap();
  _searcher->clearNodeMap();
  return true;
}

void visBSplineTrajectory(ewok::UniformBSpline3D<6, double>::Ptr &b_spline_, string ns) {
  visualization_msgs::Marker traj_vis;
  double dt = b_spline_->dt();
  traj_vis.header.stamp = ros::Time::now();
  traj_vis.header.frame_id = "world";

  traj_vis.ns = ns;
  traj_vis.id = 0;
  traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;

  traj_vis.action = visualization_msgs::Marker::DELETE;

  traj_vis.action = visualization_msgs::Marker::ADD;
  traj_vis.scale.x = _vis_traj_width;
  traj_vis.scale.y = _vis_traj_width;
  traj_vis.scale.z = _vis_traj_width;
  traj_vis.pose.orientation.x = 0.0;
  traj_vis.pose.orientation.y = 0.0;
  traj_vis.pose.orientation.z = 0.0;
  traj_vis.pose.orientation.w = 1.0;

  double traj_len = 0.0;
  int count = 0;
  Vector3d cur, pre;
  cur.setZero();
  pre.setZero();

  traj_vis.points.clear();

  Vector3d state;
  double velocity;
  geometry_msgs::Point pt;
  std_msgs::ColorRGBA colors;

  int cureve_seg_num = b_spline_->size() - b_spline_->N + 1;
  ROS_INFO("control_point num : %d ", b_spline_->size());
  for (double t = 0.0; t < cureve_seg_num * dt; t += 0.06 * dt, count += 1) {
    state = b_spline_->evaluate(t, 0);
    velocity = ((b_spline_->evaluate(t, 1)).norm() - 0.2) / 0.9;
    cur(0) = pt.x = state(0);
    cur(1) = pt.y = state(1);
    cur(2) = pt.z = state(2);

//        if(velocity<=0.3333){
//            colors.r=0.0;
//            colors.g=velocity/0.3333;
//            colors.b=1.0-velocity/0.3333;
//        }
//        else if(velocity<=0.6666){
//            colors.r=(velocity-0.3333)/0.3333;
//            colors.g=1.0;
//            colors.b=0.0;
//        }
//        else{
//            colors.r=1.0;
//            colors.g=1.0-(velocity-0.6666)/0.3333;
//            colors.b=0.0;
//        }
    colors.r = 1.0;
    colors.g = 0.0;
    colors.b = 0.0;
    colors.a = 0.6;

    traj_vis.points.push_back(pt);
    traj_vis.colors.push_back(colors);
    if (count) traj_len += (pre - cur).norm();
    pre = cur;
  }

  ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
  _traj_vis_pub.publish(traj_vis);

}
void visBSplineControlPts(ewok::UniformBSpline3D<6, double>::Ptr &b_spline_, string ns) {
  visualization_msgs::MarkerArray grid_vis;
  for (auto &mk: grid_vis.markers)
    mk.action = visualization_msgs::Marker::DELETE;

  grid_vis.markers.clear();

  visualization_msgs::Marker mk;
  mk.header.frame_id = "world";
  mk.header.stamp = ros::Time::now();
  mk.ns = ns;
  mk.type = visualization_msgs::Marker::CUBE;
  mk.action = visualization_msgs::Marker::ADD;

  mk.pose.orientation.x = 0.0;
  mk.pose.orientation.y = 0.0;
  mk.pose.orientation.z = 0.0;
  mk.pose.orientation.w = 1.0;
  mk.color.a = 1.0;
  mk.color.r = 1.0;
  mk.color.g = 0.0;
  mk.color.b = 0.0;

  int idx = 0;
  for (int i = 0; i < int(b_spline_->size()); i++) {
    mk.id = idx;
    auto ctrpt = b_spline_->getControlPoint(i);
    mk.pose.position.x = ctrpt(0);
    mk.pose.position.y = ctrpt(1);
    mk.pose.position.z = ctrpt(2);

    mk.scale.x = _resolution;
    mk.scale.y = _resolution;
    mk.scale.z = _resolution;

    idx++;
    grid_vis.markers.push_back(mk);
  }

  _grid_path_vis_pub.publish(grid_vis);
}
void visExpNode(vector<BSplineGridNodePtr> nodes) {
  visualization_msgs::Marker node_vis;
  node_vis.header.frame_id = "world";
  node_vis.header.stamp = ros::Time::now();
  node_vis.ns = "kinodynamic/visited_nodes";
  node_vis.type = visualization_msgs::Marker::CUBE_LIST;
  node_vis.action = visualization_msgs::Marker::ADD;
  node_vis.id = 0;

  node_vis.pose.orientation.x = 0.0;
  node_vis.pose.orientation.y = 0.0;
  node_vis.pose.orientation.z = 0.0;
  node_vis.pose.orientation.w = 1.0;
  node_vis.color.a = 0.3;
  node_vis.color.r = 0.0;
  node_vis.color.g = 1.0;
  node_vis.color.b = 0.0;

  node_vis.scale.x = _resolution;
  node_vis.scale.y = _resolution;
  node_vis.scale.z = _resolution;

  geometry_msgs::Point pt;
  for (int i = 0; i < int(nodes.size()); i++) {
    Vector3d coord = nodes[i]->coord;
    pt.x = coord(0);
    pt.y = coord(1);
    pt.z = coord(2);
    node_vis.points.push_back(pt);
  }

  _nodes_vis_pub.publish(node_vis);
}

void visVel() {
  visualization_msgs::Marker node_vis;
  node_vis.header.frame_id = "world";
  node_vis.header.stamp = ros::Time::now();
  node_vis.ns = "start";
  node_vis.type = visualization_msgs::Marker::CUBE_LIST;
  node_vis.action = visualization_msgs::Marker::ADD;
  node_vis.id = 0;

  node_vis.pose.orientation.x = 0.0;
  node_vis.pose.orientation.y = 0.0;
  node_vis.pose.orientation.z = 0.0;
  node_vis.pose.orientation.w = 1.0;

  node_vis.scale.x = _resolution;
  node_vis.scale.y = _resolution;
  node_vis.scale.z = _resolution;

  geometry_msgs::Point pt;
  std_msgs::ColorRGBA cl;

  pt.x = _end_pt(0);
  pt.y = _end_pt(1);
  pt.z = _end_pt(2);
  node_vis.points.push_back(pt);
  cl.a = 0.7;
  cl.r = 1.0;
  cl.g = 1.0;
  cl.b = 0.0;
  node_vis.colors.push_back(cl);

  pt.x = _end_pt(0) - _end_vel(0);
  pt.y = _end_pt(1) - _end_vel(1);
  pt.z = _end_pt(2) - _end_vel(2);
  node_vis.points.push_back(pt);
  cl.a = 0.4;
  cl.r = 1.0;
  cl.g = 1.0;
  cl.b = 0.0;
  node_vis.colors.push_back(cl);

  pt.x = _start_pt(0);
  pt.y = _start_pt(1);
  pt.z = _start_pt(2);
  node_vis.points.push_back(pt);
  cl.a = 0.7;
  cl.r = 1.0;
  cl.g = 0.0;
  cl.b = 0.0;
  node_vis.colors.push_back(cl);

  pt.x = _start_pt(0) + _start_vel(0);
  pt.y = _start_pt(1) + _start_vel(1);
  pt.z = _start_pt(2) + _start_vel(2);
  node_vis.points.push_back(pt);
  cl.a = 0.7;
  cl.r = 1.0;
  cl.g = 0.0;
  cl.b = 0.0;
  node_vis.colors.push_back(cl);

  _end_vis_pub.publish(node_vis);
}

quadrotor_msgs::BSplineTrajectory getBsplineTraj(ewok::UniformBSpline3D<6, double>::Ptr bspline) {
  quadrotor_msgs::BSplineTrajectory traj;
  traj.action = quadrotor_msgs::BSplineTrajectory::ACTION_ADD;
  traj.order = 6;
  traj.delt_t = _delt_t;
  for (int i = 0; i < int(bspline->size()); i++) {
    auto ctrpt = bspline->getControlPoint(i);
    traj.control_pt_x.push_back(ctrpt(0));
    traj.control_pt_y.push_back(ctrpt(1));
    traj.control_pt_z.push_back(ctrpt(2));
  }
  traj.header.frame_id = "/bernstein";
  traj.header.stamp = _odom.header.stamp;
  _start_time = traj.header.stamp;
  traj.trajectory_id = _traj_id;
  return traj;
}

void bnuk_log_write_VFRBK(const std::string file_name) {
  ofstream outfile(file_name);

  outfile << "bnuk_traj_duration" << ',' << bnuk_log.traj_duration << endl;
  outfile << "bnuk_traj_lenth" << ',' << bnuk_log.traj_lenth << endl;
  outfile << "computime_ms" << ',' << bnuk_log.compute_time_mean << endl;
  outfile << "replan" << ',' << bnuk_log.replan << endl;
  outfile << "traj_lenth_vec,";
  for (auto &it : bnuk_log.traj_lenth_vec)
    outfile << "," << (it);
  outfile << endl;
  outfile << "traj_time,";
  for (auto &it : bnuk_log.traj_time)
    outfile << "," << (it);
  outfile << endl;

  outfile << "-----Opt------" << endl;
  outfile << "bnuk_traj_duration" << ',' << bnuk_opt_log.traj_duration << endl;
  outfile << "bnuk_traj_lenth" << ',' << bnuk_opt_log.traj_lenth << endl;
  outfile << "computime_ms" << ',' << bnuk_opt_log.compute_time_mean << endl;
  outfile << "replan" << ',' << bnuk_opt_log.replan << endl;
  outfile << "traj_lenth_vec,";
  for (auto &it : bnuk_opt_log.traj_lenth_vec)
    outfile << "," << (it);
  outfile << endl;
  outfile << "traj_time,";
  for (auto &it : bnuk_opt_log.traj_time)
    outfile << "," << (it);
  outfile << endl;
  outfile.close();
}
bool solve_opt() {
  ROS_INFO("BNUK_OPT start");
  b_spline_optimizer optimizer;
  vector<double> QboudRadius = _searcher->getCentersRadius().second;
  vector<Vector3d> QboudCenters = _searcher->getCentersRadius().first;

  vector<Eigen::Vector3d> ControlPts = _searcher->getpath();

  vector<double> weight;
  for (int i = 0; i < ControlPts.size(); ++i) {
    weight.emplace_back(0.1);
  }

  optimizer.setMapRange(Vector3d(-25, -25, 0), Vector3d(25, 25, 5));
  optimizer.setDynamicCons(_dynamics_max_vel, _dynamics_max_acc);
  optimizer.setEndPos(_end_state[0]);
  optimizer.setEndVel(_end_state[1]);
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  optimizer.usingQcon();
  optimizer.setMQM(2);
  optimizer.setControlPts(ControlPts);
  optimizer.setWeight(weight);
  optimizer.setQobjVal();
  optimizer.setCobjVal_SafeLimit();
  optimizer.setAVal_DynamicCons();
  optimizer.setAVal_EndCons();
  optimizer.setQconVal(QboudCenters, QboudRadius);
//  optimizer.printProblem();
  ros::Time time_opt_s = ros::Time::now();
  bool result = optimizer.solve();
  ros::Time time_opt_e = ros::Time::now();
  bnuk_opt_log.compute_time.emplace_back((time_opt_e - time_opt_s).toNSec() / 1000000.0);
  bnuk_opt_log.compute_time_mean =
      accumulate(bnuk_opt_log.compute_time.begin(), bnuk_opt_log.compute_time.end(), 0.0)
          / bnuk_opt_log.compute_time.size();
  std::cout << "opt takes "
            << (time_opt_e - time_opt_s).toNSec() / 1000000.0 << "ms.\n";

  double *opt_xx = optimizer.xx;

  vector<Eigen::Vector3d> opt_cps;

  for (int i = 0; i < ControlPts.size(); ++i) {
    opt_cps.emplace_back(opt_xx[i], opt_xx[ControlPts.size() + i], opt_xx[2 * ControlPts.size() + i]);
    //cout << opt_xx[i] << "," << opt_xx[ControlPts.size() + i] << "," << opt_xx[2*ControlPts.size() + i] << endl;
  }

  cout << endl;
  _b_spline_->ClearControlPoints();
  for (auto &ControlPt : ControlPts) {
    _b_spline_->push_back(ControlPt);
    //cout << ControlPt(0) << "," << ControlPt(1) << "," << ControlPt(2) << endl;
  }

  _b_spline_opt_->ClearControlPoints();
  for (auto &opt_cp : opt_cps) {
    _b_spline_opt_->push_back(opt_cp);
  }
  if (!result) {
    ROS_ERROR("BNUK_OPT error!!!!!!!!!!!!");
    return false;
  }

}

