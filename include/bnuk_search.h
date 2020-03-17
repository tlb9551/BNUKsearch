//
// Created by tlb on 18-9-14.
//

#ifndef PROJECT_BNUK_SEARCH_H
#define PROJECT_BNUK_SEARCH_H
#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"
#include <stdio.h>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <sdf_tools/collision_map.hpp>

#define inf 1>>30

struct BSplineGridNode;
typedef BSplineGridNode *BSplineGridNodePtr;
struct BSplineGridNode {
  int id;        // 1--> open set, -1 --> closed set
  Eigen::Vector3d coord;
  Eigen::Vector3i index;
  double steplen = 0.0;
  struct {
    int has_center;//0 not visited; 1 visited has_center ; 2 visited is center (know)
    int EDSquare;//minimum euler distance square to collision
    Eigen::Vector3i center;

  } ed2c;
  double gScore, fScore;
  BSplineGridNodePtr cameFrom;
  typename std::multimap<double, BSplineGridNodePtr>::iterator nodeMapIt;
  double occupancy = 0;

  std::vector<BSplineGridNodePtr> hisNodeList; // use a list to record nodes in its history

  BSplineGridNode(Eigen::Vector3i _index) {
    id = 0;
    index = _index;

    ed2c.EDSquare = 0;
    ed2c.has_center = 0;
    steplen = 0.0;
    gScore = inf;
    fScore = inf;
    cameFrom = NULL;
  }

  BSplineGridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord) {
    id = 0;
    index = _index;
    coord = _coord;
    ed2c.EDSquare = 0;
    ed2c.has_center = 0;
    steplen = 0.0;

    gScore = inf;
    fScore = inf;
    cameFrom = NULL;
  }

  BSplineGridNode() {

    ed2c.EDSquare = 0;
    ed2c.has_center = 0;
    steplen = 0.0;
    id = 0;
    gScore = inf;
    fScore = inf;
    cameFrom = NULL;
  }

  ~BSplineGridNode() {}
};

class bnuk_search {
 public:
  const static int _DEG = 5;
  typedef Eigen::Matrix<double, 1, _DEG + 1> VectorNT;
  typedef Eigen::Matrix<double, _DEG + 1, 1> VectorN;
  typedef Eigen::Matrix<double, _DEG + 1, _DEG + 1> MatrixN;
  typedef std::array<MatrixN, 4> MatrixNArray;

  bnuk_search(Eigen::Vector3i size, double _dt) {
    x_size = size(0);
    y_size = size(1);
    z_size = size(2);
    dt = _dt;
    double current_pow_inv_dt = 1.0;
    for (int i = 0; i < 2 * (_DEG + 1); ++i) {
      pow_inv_dt_[i] = current_pow_inv_dt;
      current_pow_inv_dt /= dt;
    }
    quadratic_cost_weight.fill(0.0);
    diff_matrix = MatrixN::Zero();
    for (int i = 0; i < 5; i++)
      diff_matrix(i, i + 1) = i + 1;
    diff_matrix = diff_matrix / dt;
    bspline_basis_matrix = computeBasisMatrix_Bspline();
    bezier_basis_matrix = computeBasisMatrix_Bezier();
    bspline2bezier_matrix =
        bezier_basis_matrix.inverse() * bspline_basis_matrix;
    bspline2bezier_diff1_matrix =
        computeBasisMatrix_Bezier(1).inverse() * diff_matrix * bspline_basis_matrix;
    bspline2bezier_diff2_matrix =
        computeBasisMatrix_Bezier(2).inverse() * diff_matrix * diff_matrix * bspline_basis_matrix;

    quadratic_cost_jacobian = computeQuadraticCostJacobian();

    start_matrix = computeStartMatrix();
    end_matrix = start_matrix;
    start_matrix_inv = start_matrix.inverse();
    end_matrix_inv = start_matrix_inv;
  };
  bnuk_search() {};
  ~bnuk_search() {

  };
  void initGridNodeMap(double _resolution);
  void linkCollsionMap(std::shared_ptr<sdf_tools::CollisionMapGrid> collision_map);
  void resetLocalMap();
  void clearNodeMap();
  void V_F_RBK_search();
  void resetCollisionMap();
  void resetPath();
  std::vector<Eigen::Vector3d> getpath();
  std::pair<std::vector<Eigen::Vector3d>, std::vector<double>> getCentersRadius();
  std::vector<BSplineGridNodePtr> getVisitedNodes();
  double getSteplen(int x, int y, int z);
  double getSteplen_opt(int x, int y, int z);
  double getSteplen_edt(int x, int y, int z);
  double getMaxVel() { return max_vel; }
  void setMaxVel(double vel);
  Eigen::Vector4f getDynamicsConstraint();
  void setDynamicsConstraint(Eigen::Vector4f dynamic_constraint);
  void setStartState(std::vector<Eigen::Vector3d> state);
  void setEndState(std::vector<Eigen::Vector3d> state);
  void setTieBreaker(double t) { tie_breaker = t; };
  void setHeuCoeff(double t) { heu_coeff = t; };
  inline void setQuadraticCostWeight(std::array<double, _DEG - 1> weight_array) {
    quadratic_cost_weight = weight_array;
  };
  std::array<double, 2 * (_DEG + 1)> pow_inv_dt_;
  std::shared_ptr<sdf_tools::CollisionMapGrid> CollisionMap;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  MatrixN bspline_basis_matrix;
  MatrixN bezier_basis_matrix;
  MatrixN bspline2bezier_matrix;
  MatrixN bspline2bezier_diff1_matrix;
  MatrixN bspline2bezier_diff2_matrix;
  MatrixN diff_matrix;
  MatrixNArray quadratic_cost_jacobian;
  std::array<double, _DEG - 1> quadratic_cost_weight;

  Eigen::Matrix<double, _DEG, _DEG> start_matrix_inv;
  Eigen::Matrix<double, _DEG, _DEG> end_matrix_inv;
  Eigen::Matrix<double, _DEG, _DEG> start_matrix;
  Eigen::Matrix<double, _DEG, _DEG> end_matrix;

  Eigen::Matrix<double, _DEG, 3> start_state;
  Eigen::Matrix<double, 2, 3> end_state;

  std::vector<std::vector<Eigen::Vector3i> > distance_queue;
  void initDistanceQueue(int maxDis);
  inline int distanceSquare_2c(int px, int py, int pz){ return px * px + py * py + pz * pz;};

 //private:
  double resolution, inv_resolution;
  double tie_breaker = 1.0 + 1.0 / 10.0;
  int deg = _DEG;
  double dt = 1.0;
  double max_vel = 0.2;
  double heu_coeff = 1.0;
  double dis_est;
  double dynamics_max_vel = 3.0, dynamics_min_vel = 0.0, dynamics_max_acc = 3.0, dynamics_min_acc = 0.0;

  bool distance_queue_init = false, gridnode_map_init = false, max_vel_set = false;

  BSplineGridNodePtr ***BSplineGridNodeMap;
  BSplineGridNode start_node[_DEG];
  BSplineGridNode end_node[_DEG - 1];

  std::vector<BSplineGridNodePtr> expandedNodes;
  std::vector<BSplineGridNodePtr> gridPath;

  int x_size, y_size, z_size;
  Eigen::Vector3d map_origin;

  std::multimap<double, BSplineGridNodePtr> openSet;
  Eigen::Matrix<double, _DEG + 1, 3> local_control_points;

  Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
  Eigen::Vector3i coord2gridIndex(Eigen::Vector3d coord);
  BSplineGridNodePtr coord2gridNodePtr(Eigen::Vector3d pos);
  BSplineGridNodePtr index2gridNodePtr(Eigen::Vector3i index);

  double getDiagHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2);
  double getDiagDis(BSplineGridNodePtr node1, BSplineGridNodePtr node2);
  double getHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2);
  double getQuadraticCostHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2);

  std::vector<BSplineGridNodePtr> retrievePath(BSplineGridNodePtr current);

  double getgscore(BSplineGridNodePtr tentative);
  double quadraticCost(const MatrixN &quadratic_cost, BSplineGridNodePtr tentative);
  double quadraticCost(int deriv, BSplineGridNodePtr tentative);

  std::vector<Eigen::Vector3i> Neighbor(Eigen::Vector3i current_idx);

  const Eigen::Matrix<double, _DEG, _DEG> computeStartMatrix();
  const Eigen::Matrix<double, _DEG, _DEG> computeEndMatrix();
  void compute_end_node();
  void compute_start_node();

  const MatrixNArray computeQuadraticCostJacobian();
  const MatrixNArray computeQuadraticCoefficients();

  const MatrixN computeBasisMatrix_Bspline(int diff_degree = 0);
  const MatrixN computeBasisMatrix_Bezier(int diff_degree = 0);

  MatrixN computeBaseCoefficients();

  void setLocalControlPoints(BSplineGridNodePtr current);
  bool checkTrajDynamics(BSplineGridNodePtr neighbor);
  bool checkTrajCollision(BSplineGridNodePtr neighbor);

  inline bool checkPtCollision(int x, int y, int z);
  double velMapping(double d, double max_vel);
  uint64_t C_n_k(uint64_t n, uint64_t k);
};

#include "../src/bnuk_search.cpp"
#endif //PROJECT_BNUK_SEARCH_H
