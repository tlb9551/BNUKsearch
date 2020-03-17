#include "ros/ros.h"
#include "ros/console.h"
#include "nav_msgs/Odometry.h"
#include "quadrotor_msgs/BSplineTrajectory.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "geometry_msgs/PoseStamped.h"
#include "tf/tf.h"
#include "tf/transform_datatypes.h"
#include <eigen3/Eigen/Dense>
#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"
#include "uniform_bspline_3d.h"
//#include "../../../devel/include/quadrotor_msgs/PositionCommand.h"

ewok::UniformBSpline3D<6, double> _b_spline_(1.0);
const int _DIM_x = 0;
const int _DIM_y = 1;
const int _DIM_z = 2;

using namespace std;
using namespace Eigen;

int _poly_order_min, _poly_order_max;

class Bnuk_TrajectoryServer {
 private:

  // Subscribers
  ros::Subscriber _odom_sub;
  ros::Subscriber _traj_sub;

  // publishers
  ros::Publisher _cmd_pub;
  ros::Publisher _vis_cmd_pub;
  ros::Publisher _vis_vel_pub;
  ros::Publisher _vis_acc_pub;

  // configuration for trajectory
  int _traj_id = 0;
  uint32_t _traj_flag = 0;

  ros::Time _final_time = ros::TIME_MIN;
  ros::Time _start_time = ros::TIME_MAX;

  // state of the server
  enum ServerState { INIT, TRAJ, HOVER } state = INIT;;
  nav_msgs::Odometry _odom;
  quadrotor_msgs::PositionCommand _cmd;
  geometry_msgs::PoseStamped _vis_cmd;

  visualization_msgs::Marker _vis_vel, _vis_acc;
  Eigen::Matrix3d rot_yaw;

  double vel_t = 0.0;

 public:

  vector<Eigen::VectorXd>
      CList;  // Position coefficients vector, used to record all the pre-compute 'n choose k' combinatorial for the bernstein coefficients .
  vector<Eigen::VectorXd> CvList; // Velocity coefficients vector.
  vector<Eigen::VectorXd> CaList; // Acceleration coefficients vector.

  Bnuk_TrajectoryServer(ros::NodeHandle &handle) {
    _odom_sub =
        handle.subscribe("/odom/fake_odom", 50, &Bnuk_TrajectoryServer::rcvOdometryCallback, this);

    _traj_sub =
        handle.subscribe("/bnuk_node/trajectory", 2, &Bnuk_TrajectoryServer::rcvTrajectoryCallabck, this);

    _cmd_pub =
        handle.advertise<quadrotor_msgs::PositionCommand>("/position_cmd", 50);

    _vis_cmd_pub =
        handle.advertise<geometry_msgs::PoseStamped>("desired_position", 50);

    _vis_vel_pub =
        handle.advertise<visualization_msgs::Marker>("desired_velocity", 50);

    _vis_acc_pub =
        handle.advertise<visualization_msgs::Marker>("desired_acceleration", 50);

    double pos_gain[3] = {5.7, 5.7, 6.2};
    double vel_gain[3] = {3.4, 3.4, 4.0};
    setGains(pos_gain, vel_gain);
  }

  void setGains(double pos_gain[3], double vel_gain[3]) {
    _cmd.kx[_DIM_x] = pos_gain[_DIM_x];
    _cmd.kx[_DIM_y] = pos_gain[_DIM_y];
    _cmd.kx[_DIM_z] = pos_gain[_DIM_z];

    _cmd.kv[_DIM_x] = vel_gain[_DIM_x];
    _cmd.kv[_DIM_y] = vel_gain[_DIM_y];
    _cmd.kv[_DIM_z] = vel_gain[_DIM_z];
  }

  void rcvOdometryCallback(const nav_msgs::Odometry &odom) {
    if (odom.child_frame_id == "X" || odom.child_frame_id == "O") return;
    // #1. store the odometry
    _odom = odom;
    _vis_cmd.header = _odom.header;
    _vis_cmd.header.frame_id = "/world";

    if (state == INIT) {
      //ROS_WARN("[TRAJ SERVER] Pub initial pos command");
      _cmd.position = _odom.pose.pose.position;

      _cmd.header.stamp = _odom.header.stamp;
      _cmd.header.frame_id = "/world";
      _cmd.trajectory_flag = _traj_flag;

      _cmd.velocity.x = 0.0;
      _cmd.velocity.y = 0.0;
      _cmd.velocity.z = 0.0;

      _cmd.acceleration.x = 0.0;
      _cmd.acceleration.y = 0.0;
      _cmd.acceleration.z = 0.0;
      _cmd_pub.publish(_cmd);

      _vis_cmd.pose.position.x = _cmd.position.x;
      _vis_cmd.pose.position.y = _cmd.position.y;
      _vis_cmd.pose.position.z = _cmd.position.z;
      _vis_cmd_pub.publish(_vis_cmd);

      return;
    }

    // #2. try to calculate the new state
    if (state == TRAJ && ((odom.header.stamp - _start_time).toSec() > (_final_time - _start_time).toSec())) {
      state = HOVER;
      _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
    }
    // #3. try to publish command
    pubPositionCommand();
  }

  void rcvTrajectoryCallabck(const quadrotor_msgs::BSplineTrajectory &traj) {
    if (traj.action == quadrotor_msgs::BSplineTrajectory::ACTION_ADD) {
      ROS_WARN("[SERVER] Loading the trajectory.");
      if ((int) traj.trajectory_id < _traj_id) return;

      _b_spline_.ClearControlPoints();

      state = TRAJ;
      _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
      _traj_id = traj.trajectory_id;
      _final_time = _start_time = traj.header.stamp;
      for (int idx = 0; idx < traj.control_pt_x.size(); idx++) {
        _b_spline_.push_back(Vector3d(traj.control_pt_x[idx], traj.control_pt_y[idx], traj.control_pt_z[idx]));
      }
      _final_time += ros::Duration(_b_spline_.maxValidTime());

    } else if (traj.action == quadrotor_msgs::BSplineTrajectory::ACTION_ABORT) {
      ROS_WARN("[SERVER] Aborting the trajectory.");
      state = HOVER;
      _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
    } else if (traj.action == quadrotor_msgs::BSplineTrajectory::ACTION_WARN_IMPOSSIBLE) {
      state = HOVER;
      _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_IMPOSSIBLE;
    }
  }

  void pubPositionCommand() {
    if (state == INIT) return;
    if (state == HOVER) {
      if (_cmd.header.frame_id != "/world") {
        _cmd.position = _odom.pose.pose.position;
      }

      _cmd.header.stamp = _odom.header.stamp;
      _cmd.header.frame_id = "/world";
      _cmd.trajectory_flag = _traj_flag;

      _cmd.velocity.x = 0.0;
      _cmd.velocity.y = 0.0;
      _cmd.velocity.z = 0.0;

      _cmd.acceleration.x = 0.0;
      _cmd.acceleration.y = 0.0;
      _cmd.acceleration.z = 0.0;
    }

    if (state == TRAJ) {
      _cmd.header.stamp = _odom.header.stamp;

      _cmd.header.frame_id = "/world";
      _cmd.trajectory_flag = _traj_flag;
      _cmd.trajectory_id = _traj_id;

      double t = max(0.0, (_odom.header.stamp - _start_time).toSec());

      _cmd.yaw_dot = 0.0;
      _cmd.yaw = 0.0;
      Vector3d traj_pos = _b_spline_.evaluate(t, 0);
      _cmd.position.x = traj_pos(0);
      _cmd.position.y = traj_pos(1);
      _cmd.position.z = traj_pos(2);
      Vector3d traj_vel = _b_spline_.evaluate(t, 1);
      _cmd.velocity.x = traj_vel(0);
      _cmd.velocity.y = traj_vel(1);
      _cmd.velocity.z = traj_vel(2);
      Vector3d traj_acc = _b_spline_.evaluate(t, 2);
      _cmd.acceleration.x = traj_acc(0);
      _cmd.acceleration.y = traj_acc(1);
      _cmd.acceleration.z = traj_acc(2);
    }

    _cmd_pub.publish(_cmd);

    _vis_cmd.header = _cmd.header;
    _vis_cmd.pose.position.x = _cmd.position.x;
    _vis_cmd.pose.position.y = _cmd.position.y;
    _vis_cmd.pose.position.z = _cmd.position.z;

    tf::Quaternion q_ = tf::createQuaternionFromYaw(_cmd.yaw);
    geometry_msgs::Quaternion odom_quat;
    tf::quaternionTFToMsg(q_, odom_quat);
    _vis_cmd.pose.orientation = odom_quat;
    _vis_cmd_pub.publish(_vis_cmd);

    _vis_vel.ns = "vel";
    _vis_vel.id = 0;
    _vis_vel.header.frame_id = "/world";
    _vis_vel.type = visualization_msgs::Marker::ARROW;
    _vis_vel.action = visualization_msgs::Marker::ADD;
    _vis_vel.color.a = 1.0;
    _vis_vel.color.r = 1.0;
    _vis_vel.color.g = 1.0;
    _vis_vel.color.b = 0.0;

    _vis_vel.header.stamp = _odom.header.stamp;

    _vis_vel.points.clear();
    geometry_msgs::Point pt;
    pt.x = _cmd.position.x;
    pt.y = _cmd.position.y;
    pt.z = _cmd.position.z;

    _vis_vel.points.push_back(pt);

    pt.x = _cmd.position.x + _cmd.velocity.x;
    pt.y = _cmd.position.y + _cmd.velocity.y;
    pt.z = _cmd.position.z + _cmd.velocity.z;

    _vis_vel.points.push_back(pt);

    _vis_vel.scale.x = 0.1;
    _vis_vel.scale.y = 0.2;
    _vis_vel.scale.z = 0.2;

    _vis_vel_pub.publish(_vis_vel);

    _vis_acc.ns = "acc";
    _vis_acc.id = 0;
    _vis_acc.header.frame_id = "/world";
    _vis_acc.type = visualization_msgs::Marker::ARROW;
    _vis_acc.action = visualization_msgs::Marker::ADD;
    _vis_acc.color.a = 1.0;
    _vis_acc.color.r = 1.0;
    _vis_acc.color.g = 0.0;
    _vis_acc.color.b = 0.0;

    _vis_acc.header.stamp = _odom.header.stamp;

    _vis_acc.points.clear();
    pt.x = _cmd.position.x;
    pt.y = _cmd.position.y;
    pt.z = _cmd.position.z;

    _vis_acc.points.push_back(pt);

    pt.x = _cmd.position.x + _cmd.acceleration.x;
    pt.y = _cmd.position.y + _cmd.acceleration.y;
    pt.z = _cmd.position.z + _cmd.acceleration.z;

    _vis_acc.points.push_back(pt);

    _vis_acc.scale.x = 0.1;
    _vis_acc.scale.y = 0.2;
    _vis_acc.scale.z = 0.2;

    _vis_acc_pub.publish(_vis_acc);
  }
};

int main(int argc, char **argv) {
  ros::init(argc, argv, "bnuk_trajectory_server_node");
  ros::NodeHandle handle("~");

  Bnuk_TrajectoryServer server(handle);

  sleep(1);
  ros::spin();

  return 0;
}

