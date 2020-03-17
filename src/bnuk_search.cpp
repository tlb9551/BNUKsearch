//
// Created by tlb on 18-9-14.
//
#ifndef PROJECT_BNUK_SEARCH_IMPL
#define PROJECT_BNUK_SEARCH_IMPL

#include "bnuk_search.h"

using namespace std;
using namespace Eigen;
using namespace sdf_tools;

/******private******/


Vector3d bnuk_search::gridIndex2coord(Vector3i index) {
  Vector3d pt;

  pt(0) = ((double) index(0) + 0.5) * resolution + map_origin(0);
  pt(1) = ((double) index(1) + 0.5) * resolution + map_origin(1);
  pt(2) = ((double) index(2) + 0.5) * resolution + map_origin(2);

  return pt;
}


Vector3i bnuk_search::coord2gridIndex(Vector3d coord) {
  Vector3i idx;
//    idx <<  min( max( int( (coord(0) ) * inv_resolution), 0), x_size - 1),
//            min( max( int( (coord(1) ) * inv_resolution), 0), y_size - 1),
//            min( max( int( (coord(2) ) * inv_resolution), 0), z_size - 1);
  idx << int((coord(0) - map_origin(0)) * inv_resolution),
      int((coord(1) - map_origin(1)) * inv_resolution),
      int((coord(2) - map_origin(2)) * inv_resolution);

  return idx;
}


BSplineGridNodePtr bnuk_search::coord2gridNodePtr(Vector3d pos) {
  Vector3i index = bnuk_search::coord2gridIndex(pos);
  return BSplineGridNodeMap[index(0)][index(1)][index(2)];
}

inline BSplineGridNodePtr bnuk_search::index2gridNodePtr(Vector3i index) {
  return BSplineGridNodeMap[index(0)][index(1)][index(2)];
}


double bnuk_search::getDiagHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2) {
  int dx = abs(node1->index(0) - node2->index(0));
  int dy = abs(node1->index(1) - node2->index(1));
  int dz = abs(node1->index(2) - node2->index(2));

  int diag = min(min(dx, dy), dz);
  dx -= diag;
  dy -= diag;
  dz -= diag;
  double h = 0;
  double span = floor(max_vel * inv_resolution);
  if (dx == 0) {
    h = ceil(sqrt(3.0) * diag / span) + ceil(sqrt(2.0) * min(dy, dz) / span) + ceil(abs(dy - dz) / span);
  }
  if (dy == 0) {
    h = ceil(sqrt(3.0) * diag / span) + ceil(sqrt(2.0) * min(dx, dz) / span) + ceil(abs(dx - dz) / span);
  }
  if (dz == 0) {
    h = ceil(sqrt(3.0) * diag / span) + ceil(sqrt(2.0) * min(dx, dy) / span) + ceil(abs(dx - dy) / span);
  }
//    if (dx == 0) {
//        h = ceil(sqrt(3.0)*diag/span + sqrt(2.0)*min(dy, dz)/span + abs(dy - dz)/span);
//    }
//    if (dy == 0) {
//        h = ceil(sqrt(3.0) * diag/span + sqrt(2.0) * min(dx, dz)/span + abs(dx - dz)/span);
//    }
//    if (dz == 0) {
//        h = ceil(sqrt(3.0) * diag/span + sqrt(2.0) * min(dx, dy)/span + abs(dx - dy)/span);
//    }
  //h += 1.0-getSteplen_edt(node1->index(0),node1->index(1),node1->index(2))/max_vel;

  return h;
}


double bnuk_search::getDiagDis(BSplineGridNodePtr node1, BSplineGridNodePtr node2) {
  int dx = abs(node1->index(0) - node2->index(0));
  int dy = abs(node1->index(1) - node2->index(1));
  int dz = abs(node1->index(2) - node2->index(2));

  int diag = min(min(dx, dy), dz);
  dx -= diag;
  dy -= diag;
  dz -= diag;
  double d = 0;
  if (dx == 0) {
    d = sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + abs(dy - dz);
  }
  if (dy == 0) {
    d = sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + abs(dx - dz);
  }
  if (dz == 0) {
    d = sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + abs(dx - dy);
  }
  return d;
}


double bnuk_search::getHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2) {
  return tie_breaker * bnuk_search::getDiagHeu(node1, node2);
}


std::vector<BSplineGridNodePtr> bnuk_search::retrievePath(BSplineGridNodePtr current) {
  vector<BSplineGridNodePtr> path;
  path.push_back(current);
  while (current->cameFrom != NULL) {
    current = current->cameFrom;
    path.push_back(current);
  }
  return path;
}

double bnuk_search::getQuadraticCostHeu(BSplineGridNodePtr node1, BSplineGridNodePtr node2) {
  return 0;
}


double bnuk_search::getgscore(BSplineGridNodePtr tentative) {
  double cost(0.0);
  for (int deriv = 1; deriv < _DEG; deriv++) {
    cost = cost + quadratic_cost_weight[deriv - 1] * quadraticCost(deriv, tentative);
  }
  return cost;
}


double bnuk_search::quadraticCost(int deriv, BSplineGridNodePtr tentative) {
  return quadraticCost(quadratic_cost_jacobian[deriv - 1] * pow_inv_dt_[2 * deriv - 1], tentative);
}

double bnuk_search::quadraticCost(const bnuk_search::MatrixN &quadratic_cost,
                                        BSplineGridNodePtr tentative) {
  int controlpts_num = 1;
  vector<Vector3d> controlpts_coord;
  BSplineGridNodePtr it = tentative;
  for (controlpts_coord.push_back(it->coord); it->cameFrom != NULL;) {
    controlpts_coord.push_back(it->coord);
    controlpts_num++;
    it = it->cameFrom;
    if (controlpts_num == _DEG + 1) break;
  }
  if (controlpts_coord.size() < _DEG + 1) return 0;

  double inc_cost[3];
  for (int i = 0; i < 3; i++) {
    VectorN coefficients;
    for (int j = 0; j < _DEG + 1; j++)
      coefficients[j] = controlpts_coord[j](i);
    inc_cost[i] =
        coefficients.transpose() * quadratic_cost
            * coefficients;
  }
  return inc_cost[0] + inc_cost[1] + inc_cost[2] + tentative->gScore;

}

const Matrix<double, 5, 5> bnuk_search::computeStartMatrix() {
  Matrix<double, _DEG, _DEG> startm;
  for (int i = 0; i < _DEG; i++) {
    for (int j = 0; j < _DEG; j++) {
      if (i == j) {
        int factorial = 1;
        for (int k = j; k > 0; k--)
          factorial *= k;
        startm(i, j) = factorial * 1.0 / pow_inv_dt_[j];
      } else startm(i, j) = 0.0;
    }
  }
  return (startm * bspline_basis_matrix.block(0, 0, _DEG, _DEG));
}

const Matrix<double, 5, 5> bnuk_search::computeEndMatrix() {
  Matrix<double, _DEG + 1, _DEG + 1> endm;
  for (int i = 0; i < _DEG + 1; i++) {
    for (int j = 0; j < _DEG + 1; j++) {
      if (j >= i) {
        int factorial = 1;
        for (int k = 0; k < i; k++)
          factorial *= (j - k);
        endm(i, j) = factorial * 1.0 / pow_inv_dt_[j];
      } else endm(i, j) = 0.0;
    }
  }

  endm = endm * bspline_basis_matrix;
  return endm.block(0, 1, _DEG, _DEG);
}

void bnuk_search::compute_end_node() {
  Matrix<double, 2, 3> p_56;
  Matrix<double, 2, 3> PV = end_state;
  Matrix<double, 3, 3> p_234;

  p_234.block(0, 0, 1, 3) = end_node[0].cameFrom->coord.transpose();
  p_234.block(1, 0, 1, 3) = end_node[0].coord.transpose();
  p_234.block(2, 0, 1, 3) = end_node[1].coord.transpose();

  p_56 = end_matrix.block(0, 3, 2, 2).inverse() * (PV - end_matrix.block(0, 0, 2, 3) * p_234);

  end_node[2].coord = p_56.block(0, 0, 1, 3).transpose();
  end_node[2].index = coord2gridIndex(end_node[2].coord);
  end_node[2].cameFrom = &end_node[1];

  end_node[3].coord = p_56.block(1, 0, 1, 3).transpose();
  end_node[3].index = coord2gridIndex(end_node[3].coord);
  end_node[3].cameFrom = &end_node[2];

}


const typename bnuk_search::MatrixN bnuk_search::computeBasisMatrix_Bspline(int diff_degree) {
  if (diff_degree > _DEG) {
    cerr << "bnuk_search::MatrixN bnuk_search::computeBasisMatrix_Bspline: diff_degree  = " << diff_degree
         << endl;
    return bnuk_search::MatrixN::Zero();
  }

  typename bnuk_search::MatrixN m;
  m.setZero();
  for (int i = 0; i < _DEG + 1 - diff_degree; ++i) {
    for (int j = 0; j < _DEG + 1 - diff_degree; ++j) {
      double sum = 0;

      for (int s = j; s < _DEG + 1 - diff_degree; ++s) {
        sum += std::pow(-1.0, s - j) * C_n_k(_DEG + 1 - diff_degree, s - j)
            * std::pow(_DEG + 1 - diff_degree - s - 1.0, _DEG + 1 - diff_degree - 1.0 - i);
      }
      m(i, j) = C_n_k(_DEG + 1 - diff_degree - 1, _DEG + 1 - diff_degree - 1 - i) * sum;
    }
  }

  uint64_t factorial = 1.0;
  for (int i = 2; i < _DEG + 1 - diff_degree; ++i) {
    factorial *= i;
  }
  m = m / factorial;
  if (diff_degree > 0) {
    for (int i = 0; i < diff_degree; i++)
      m(_DEG - i, _DEG - i) = 1.0;
  }
  return m.template cast<double>();

}


const typename bnuk_search::MatrixN bnuk_search::computeBasisMatrix_Bezier(int diff_degree) {
  if (diff_degree > _DEG) {
    cerr << "bnuk_search::MatrixN bnuk_search::computeBasisMatrix_Bezier: diff_degree  = " << diff_degree
         << endl;
    return bnuk_search::MatrixN::Zero();
  }

  typename bnuk_search::MatrixN m;
  m.setZero();

  for (int i = 0; i < _DEG + 1 - diff_degree; ++i) {
    for (int j = 0; j < _DEG + 1 - diff_degree; ++j) {
      if (i < j)
        m(i, j) = 0;
      else {
        m(i, j) = std::pow(-1.0, i - j) * C_n_k(_DEG - diff_degree, j) * C_n_k(_DEG - diff_degree - j, i - j);
      }
    }
  }

  if (diff_degree > 0) {
    for (int i = 0; i < diff_degree; i++)
      m(_DEG - i, _DEG - i) = 1.0;
  }

  return (m).template cast<double>();
}

const typename bnuk_search::MatrixNArray bnuk_search::computeQuadraticCostJacobian() {
  typename bnuk_search::MatrixNArray m_array;
  typename bnuk_search::MatrixN bspline_basis_matrix;
  bspline_basis_matrix = bnuk_search::computeBasisMatrix_Bspline();
  typename bnuk_search::MatrixNArray
      quadratic_coefficients = bnuk_search::computeQuadraticCoefficients();

  for (int i = 0; i < _DEG - 1; ++i) {
    bnuk_search::MatrixN m;
    m = bspline_basis_matrix.transpose() *
        quadratic_coefficients[i] * bspline_basis_matrix;

    m_array[i] = m.template cast<double>();
  }
  return m_array;
}

const typename bnuk_search::MatrixNArray bnuk_search::computeQuadraticCoefficients() {
  typename bnuk_search::MatrixNArray res;
  typename bnuk_search::MatrixN base_coefficients =
      bnuk_search::computeBaseCoefficients();

  for (int derivative = 1; derivative < _DEG; derivative++) {
    res[derivative - 1].setZero();

    for (int col = 0; col < _DEG + 1 - derivative; col++) {
      for (int row = 0; row < _DEG + 1 - derivative; row++) {
        double exp = (_DEG - derivative) * 2 + 1 - row - col;

        res[derivative - 1](_DEG - row, _DEG - col) =
            base_coefficients(derivative, _DEG + 1 - 1 - row) *
                base_coefficients(derivative, _DEG + 1 - 1 - col) * 2 / exp;
      }
    }
  }

  return res;
}

typename bnuk_search::MatrixN bnuk_search::computeBaseCoefficients() {
  typename bnuk_search::MatrixN base_coefficients;

  base_coefficients.setZero();
  base_coefficients.row(0).setOnes();

  const int DEG = _DEG;
  int order = DEG;
  for (int n = 1; n < _DEG + 1; n++) {
    for (int i = DEG - order; i < _DEG + 1; i++) {
      base_coefficients(n, i) = (order - DEG + i) * base_coefficients(n - 1, i);
    }
    order--;
  }
  return base_coefficients.template cast<double>();
}

void bnuk_search::setLocalControlPoints(BSplineGridNodePtr current) {
  for (int i = 0; i < _DEG; i++) {
    local_control_points.row(_DEG - 1 - i) = current->coord.transpose();
    if (current->cameFrom != NULL)
      current = current->cameFrom;
    else {
      if (i != _DEG - 1)
        cerr << "error: bnuk_search::setLocalControlPoints : i=" << i << endl;
    }
  }
}


bool bnuk_search::checkTrajDynamics(BSplineGridNodePtr neighbor) {

  if (_DEG < 2) {
    cerr << "bnuk_search::checkTrajDynamics : degree of bspline is too low" << endl;
    return false;
  }
  Matrix<double, _DEG + 1, 3> bspline_diff0_control_points;
  Matrix<double, _DEG, 3> bezier_diff1_control_points;
  Matrix<double, _DEG - 1, 3> bezier_diff2_control_points;

  bspline_diff0_control_points = local_control_points;
  bspline_diff0_control_points.row(_DEG) = neighbor->coord.transpose();

  bezier_diff1_control_points = (bspline2bezier_diff1_matrix * bspline_diff0_control_points).block(0, 0, _DEG, 3);

  bezier_diff2_control_points = (bspline2bezier_diff2_matrix * bspline_diff0_control_points).block(0, 0, _DEG - 1, 3);

  if (bezier_diff1_control_points.maxCoeff() > dynamics_max_vel ||
      bezier_diff1_control_points.minCoeff() < dynamics_min_vel ||
      bezier_diff2_control_points.maxCoeff() > dynamics_max_acc ||
      bezier_diff2_control_points.minCoeff() < dynamics_min_acc
      ) {
    return false;//Not dynamically feasible
  }

  return true;//Dynamically feasible
}

bool bnuk_search::checkTrajCollision(BSplineGridNodePtr neighbor) {
  Matrix<double, _DEG + 1, 3> bspline_control_points;
  Matrix<double, _DEG + 1, 3> bezier_control_points;

  bspline_control_points = local_control_points;
  bspline_control_points.row(_DEG) = neighbor->coord.transpose();
  bezier_control_points = bspline2bezier_matrix * bspline_control_points;

  return true;
}


double bnuk_search::velMapping(double d, double max_vel) {
  return std::max(std::min(d, max_vel), 0.0);
}

std::vector<Vector3i> bnuk_search::Neighbor(Vector3i current_idx) {
  std::vector<Vector3i> nbrs;
  Vector3i nbr(0, 0, 0);

  //double current_vel=index2gridNodePtr(current_idx)->steplen;
  double current_vel = getSteplen(current_idx(0), current_idx(1), current_idx(2));

  int dir_0_0 = floor(current_vel);
  int dir_45_0 = floor(current_vel / 1.414);
  int dir_45_45 = floor(current_vel / 1.732);
//     dir_45_0=floor(current_vel);
//     dir_45_45=floor(current_vel);


  for (int i = -1; i < 2; ++i) {
    for (int j = -1; j < 2; ++j) {
      for (int k = -1; k < 2; ++k) {
        switch (abs(i) + abs(j) + abs(k)) {
          case 0:continue;
          case 1:nbr(0) = current_idx(0) + i * dir_0_0;
            nbr(1) = current_idx(1) + j * dir_0_0;
            nbr(2) = current_idx(2) + k * dir_0_0;
            break;
          case 2:nbr(0) = current_idx(0) + i * dir_45_0;
            nbr(1) = current_idx(1) + j * dir_45_0;
            nbr(2) = current_idx(2) + k * dir_45_0;
            break;
          case 3:nbr(0) = current_idx(0) + i * dir_45_45;
            nbr(1) = current_idx(1) + j * dir_45_45;
            nbr(2) = current_idx(2) + k * dir_45_45;
            break;
          default:continue;
        }
        nbrs.push_back(nbr);
      }
    }
  }
  return nbrs;
}

void bnuk_search::setStartState(vector<Vector3d> state) {
  /**** set startstate ****/
  if (state.size() > _DEG) {
    cerr << "V_F_RBK_search: startstate too much" << endl;
    return;
  }
  if (state.size() < _DEG - 1)
    cerr << "V_F_RBK_search: startstate too few" << endl;

  for (int i = 0; i < state.size(); i++) {
    start_state(i, 0) = state[i](0);
    start_state(i, 1) = state[i](1);
    start_state(i, 2) = state[i](2);
  }
  for (int i = state.size(); i < _DEG; i++) {
    start_state(i, 0) = 0;
    start_state(i, 1) = 0;
    start_state(i, 2) = 0;
  }
  /***** set start grid ****/
  Matrix<double, _DEG, 3> start_node_coord = start_matrix_inv * start_state;
  start_node[0].coord = (start_node_coord.block(0, 0, 1, 3)).transpose();
  start_node[0].index = coord2gridIndex(start_node[0].coord);
  start_node[0].cameFrom = NULL;
  for (int i = 1; i < _DEG; i++) {
    start_node[i].coord = (start_node_coord.block(i, 0, 1, 3)).transpose();
    start_node[i].index = coord2gridIndex(start_node[i].coord);
    start_node[i].cameFrom = &start_node[i - 1];
    start_node[i].gScore = 0;
  }
}

void bnuk_search::setEndState(vector<Vector3d> state) {
  if (!max_vel_set) {
    cerr << "setEndState(): set maxvel first!" << endl;
  }

  if (state.size() != 2)
    cerr << "V_F_RBK_search: endstate :" << state.size() << endl;

  for (int i = 0; i < state.size(); i++) {
    end_state(i, 0) = state[i](0);
    end_state(i, 1) = state[i](1);
    end_state(i, 2) = state[i](2);
  }

  /***** set end grid ****/
  end_node[0].coord = (end_state.block(0, 0, 1, 3) - end_state.block(1, 0, 1, 3)).transpose();
  end_node[0].index = coord2gridIndex(end_node[0].coord);
  end_node[0].gScore = 0;

  end_node[1].coord = end_state.block(0, 0, 1, 3).transpose();
  end_node[1].index = coord2gridIndex(end_node[1].coord);
  end_node[1].gScore = 0;
  end_node[1].cameFrom = &end_node[0];


  /***** set end velocity field ****/

  for (int x = -int(max_vel * inv_resolution); x <= int(max_vel * inv_resolution); x++) {
    for (int y = -int(max_vel * inv_resolution); y <= int(max_vel * inv_resolution); y++) {
      for (int z = -int(max_vel * inv_resolution); z <= int(max_vel * inv_resolution); z++) {
        if (end_node[0].index(0) + x >= 0 && end_node[0].index(0) + x < x_size
            && end_node[0].index(1) + y >= 0 && end_node[0].index(1) + y < y_size
            && end_node[0].index(2) + z >= 0 && end_node[0].index(2) + z < z_size) {
          BSplineGridNodeMap[end_node[0].index(0) + x][end_node[0].index(1) + y][end_node[0].index(2)
              + z]->ed2c.has_center = 2;
          BSplineGridNodeMap[end_node[0].index(0) + x][end_node[0].index(1) + y][end_node[0].index(2)
              + z]->ed2c.EDSquare = distanceSquare_2c(x, y, z) + 0.05;
        }
      }
    }
  }
}
/******public********/

void bnuk_search::initGridNodeMap(double _resolution) {
  gridnode_map_init = true;

  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;

  BSplineGridNodeMap = new BSplineGridNodePtr **[x_size];
  for (int i = 0; i < x_size; i++) {
    BSplineGridNodeMap[i] = new BSplineGridNodePtr *[y_size];
    for (int j = 0; j < y_size; j++) {
      BSplineGridNodeMap[i][j] = new BSplineGridNodePtr[z_size];
      for (int k = 0; k < z_size; k++) {
        Vector3i tmpIdx(i, j, k);
        BSplineGridNodeMap[i][j][k] = new BSplineGridNode(tmpIdx);
      }
    }
  }
}


void bnuk_search::linkCollsionMap(std::shared_ptr<sdf_tools::CollisionMapGrid> collision_map) {
  CollisionMap = collision_map;
  map_origin = collision_map->GetOriginTransform().translation();
  if (x_size != collision_map->GetNumXCells() ||
      y_size != collision_map->GetNumYCells() ||
      z_size != collision_map->GetNumZCells()) {
    std::cerr << "error in linkCollisionMap: size_error" << endl;
    return;
  }
//    ros::Time edt1=ros::Time::now();
//    auto EDT = collision_map->ExtractDistanceField(INFINITY);
  ros::Time edt2 = ros::Time::now();
//    ROS_INFO("ETD time cost : %f",(edt2-edt1).toNSec()/1000000.0);
  Vector3d coord;
  for (int64_t i = 0; i < x_size; i++) {
    for (int64_t j = 0; j < y_size; j++) {
      for (int64_t k = 0; k < z_size; k++) {
        coord(0) = map_origin(0) + (double) (i + 0.5) * resolution;
        coord(1) = map_origin(1) + (double) (j + 0.5) * resolution;
        coord(2) = map_origin(2) + (double) (k + 0.5) * resolution;
        Vector3i index = coord2gridIndex(coord);

        if (index(0) >= x_size || index(1) >= y_size || index(2) >= z_size
            || index(0) < 0 || index(1) < 0 || index(2) < 0)
          continue;

        BSplineGridNodePtr ptr = index2gridNodePtr(Vector3i(i, j, k));
        ptr->id = 0;
        ptr->occupancy = collision_map->Get(i, j, k).first.occupancy;
        ptr->coord = coord;
//                BSplineGridNodeMap[i][j][k]->steplen=velMapping(sqrt(EDT.GetImmutable(index).first.distance_square),max_vel*inv_resolution);

      }//index:w.r.t global_grid_map   i,j,k: w.r.t local_grid_map
    }
  }
  ros::Time map2 = ros::Time::now();
  ROS_INFO("Nodemap   time cost : %f", (map2 - edt2).toNSec() / 1000000.0);

}

void bnuk_search::resetLocalMap() {
  //ROS_WARN("expandedNodes size : %d", expandedNodes.size());
  for (auto tmpPtr:expandedNodes) {
    tmpPtr->occupancy = 0; // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  for (auto ptr:openSet) {
    BSplineGridNodePtr tmpPtr = ptr.second;
    tmpPtr->occupancy = 0; // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  expandedNodes.clear();
  //ROS_WARN("local map reset finish");
}

void bnuk_search::clearNodeMap() {
  for (int i = 0; i < x_size; i++) {
    for (int j = 0; j < y_size; j++) {
      for (int k = 0; k < z_size; k++) {
        delete BSplineGridNodeMap[i][j][k];
      }
      delete[] BSplineGridNodeMap[i][j];
    }
    delete[] BSplineGridNodeMap[i];
  }
  delete[] BSplineGridNodeMap;
}


void bnuk_search::V_F_RBK_search() {
  ros::Time time_1 = ros::Time::now();
  openSet.clear();
  BSplineGridNodePtr neighborPtr = NULL;
  BSplineGridNodePtr current = NULL;
  BSplineGridNodePtr startPtr = &start_node[_DEG - 1];
  BSplineGridNodePtr endPtr = &end_node[0];
  vector<Vector3i> neighbors_index;

  dis_est = getDiagDis(startPtr, endPtr);
  startPtr->gScore = 0;
  startPtr->fScore = getHeu(startPtr, endPtr);
  startPtr->id = 1;
  openSet.insert(make_pair(startPtr->fScore, startPtr));

  double tentative_gScore;

  int num_iter = 0;
  while (!openSet.empty()) {
    num_iter++;
    current = openSet.begin()->second;
    if (current->index == endPtr->index) {
      ROS_INFO("VF_RBK_Search reach goal...");
      ROS_INFO("total num of iteration: %d ", num_iter);
      ros::Time time_2 = ros::Time::now();
      ROS_INFO("VF_RBK_Search total time cost: %f", (time_2 - time_1).toNSec() / 1000000.0);

//    gridPath=retrievePath(current);
      end_node[0].cameFrom = current->cameFrom;

      compute_end_node();
      gridPath = retrievePath(&end_node[3]);

      return;
    }
    openSet.erase(openSet.begin());
    current->id = -1;
    expandedNodes.push_back(current);

    setLocalControlPoints(current);
    neighbors_index = Neighbor(current->index);
    for (int i = 0; i < neighbors_index.size(); i++) {

      if (neighbors_index[i](0) < 0 || neighbors_index[i](0) >= x_size
          || neighbors_index[i](1) < 0 || neighbors_index[i](1) >= y_size
          || neighbors_index[i](2) < 0 || neighbors_index[i](2) >= z_size)
        continue;
      neighborPtr = index2gridNodePtr(neighbors_index[i]);
      if (neighborPtr->occupancy > 0.5)
        continue;
      if (neighborPtr->id == -1)
        continue;
      if (!checkTrajDynamics(neighborPtr)) {
        // cout<<"NOT Dynamic feasibility "<<endl;
        continue;
      }


//            double static_cost=sqrt(distanceSquare_2c(neighbors_index[i](0)-current->index(0),
//                                                      neighbors_index[i](1)-current->index(1),
//                                                      neighbors_index[i](2)-current->index(2)));
      double static_cost = 1.0;
      tentative_gScore = current->gScore + static_cost + getgscore(neighborPtr);
      if (neighborPtr->id != 1) {
        //discover a new node
        neighborPtr->id = 1;
        neighborPtr->cameFrom = current;
        neighborPtr->gScore = tentative_gScore;
        neighborPtr->fScore = neighborPtr->gScore + getHeu(neighborPtr, endPtr);
        neighborPtr->nodeMapIt =
            openSet.insert(make_pair(neighborPtr->fScore, neighborPtr)); //put neighbor in open set and record it.

        //std::cout<<"gscore"<<tentative_gScore<<"heu"<<getHeu(neighborPtr, endPtr)<<endl;

        continue;
      } else if (tentative_gScore <= neighborPtr->gScore) {
        //in open set and need update
        neighborPtr->cameFrom = current;
        neighborPtr->gScore = tentative_gScore;
        neighborPtr->fScore = tentative_gScore + getHeu(neighborPtr, endPtr);
        openSet.erase(neighborPtr->nodeMapIt);
        neighborPtr->nodeMapIt =
            openSet.insert(make_pair(neighborPtr->fScore, neighborPtr)); //put neighbor in open set and record it.
      }
    }

  }
  ros::Time time_2 = ros::Time::now();
  ROS_WARN("Time consume in VF_RBK_searching  is %f", (time_2 - time_1).toSec());
}


void bnuk_search::resetCollisionMap() {
  //ROS_WARN("expandedNodes size : %d", expandedNodes.size());
  for (auto tmpPtr:expandedNodes) {
    tmpPtr->occupancy = 0; // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  for (auto ptr:openSet) {
    auto tmpPtr = ptr.second;
    tmpPtr->occupancy = 0; // forget the occupancy
    tmpPtr->id = 0;
    tmpPtr->cameFrom = NULL;
    tmpPtr->gScore = inf;
    tmpPtr->fScore = inf;
  }

  expandedNodes.clear();
  //ROS_WARN("local map reset finish");
}


void bnuk_search::resetPath() {
  gridPath.clear();
}


std::vector<Vector3d> bnuk_search::getpath() {
  vector<Vector3d> path;

  for (auto ptr:gridPath) {
    path.push_back(ptr->coord);
  }

  reverse(path.begin(), path.end());
  return path;
}


std::pair<std::vector<Eigen::Vector3d>, std::vector<double>> bnuk_search::getCentersRadius() {
  vector<Vector3d> centers;
  vector<double> radius;
  vector<Vector3d> control_pts;
  for(int i = gridPath.size()-1;i>=0;i--){
    control_pts.push_back(gridPath[i]->coord);
  }

  for (int i = 2;i <= control_pts.size()-4;i++) {
    Vector3d center = (control_pts[i] + control_pts[i+1])/2.0;
    Vector3d gap = control_pts[i] - control_pts[i+1];
    double min_radius = gap.norm();
    auto index = coord2gridIndex(center);

    centers.push_back(center);

    BSplineGridNodeMap[index(0)][index(1)][index(2)]->ed2c.has_center = 0;
    double ten_radius = resolution * getSteplen(index(0), index(1), index(2));
    if(ten_radius > min_radius){
      radius.push_back(max(0.0,ten_radius -1.732 * resolution / 2.0));
    }else{
      radius.push_back(max(0.0,min_radius -1.732 * resolution / 2.0));
    }
  }
  return {centers, radius};
}


std::vector<BSplineGridNodePtr> bnuk_search::getVisitedNodes() {
  vector<BSplineGridNodePtr> visited_nodes;
  for (int i = 0; i < x_size; i++)
    for (int j = 0; j < y_size; j++)
      for (int k = 0; k < z_size; k++) {
        if (BSplineGridNodeMap[i][j][k]->id != 0)
          //if(GridNodeMap[i][j][k]->id == -1)
          visited_nodes.push_back(BSplineGridNodeMap[i][j][k]);
      }
  ROS_WARN("visited_nodes size : %d", visited_nodes.size());
  return visited_nodes;
}

double bnuk_search::getSteplen(int x, int y, int z) {
  assert(!(x < 0 || x >= x_size || y < 0 || y >= y_size || z < 0 || z >= z_size));

  bool collision_detected = false;

  int ptx, pty, ptz;
  int edsquare = 0;
  for (int ds = 0; ds < distance_queue.size(); ds++) {
    if (!distance_queue[ds].empty()) {
      for (auto &it:distance_queue[ds]) {
        ptx = x + (it)(0);
        pty = y + (it)(1);
        ptz = z + (it)(2);

        if (ptx < 0 || ptx >= x_size || pty < 0 || pty >= y_size || ptz < 0 || ptz >= z_size) {//detected bounding
          collision_detected = true;
          edsquare = ds;
          break;
        }
        if (checkPtCollision(ptx, pty, ptz)) {
          collision_detected = true;
          edsquare = ds;
          break;
        }
      }
    }
    if (collision_detected)
      break;
  }
  if (!collision_detected) {
    edsquare = distance_queue.size() - 1;
  }

  if (BSplineGridNodeMap[x][y][z]->ed2c.has_center == 2){//2 visited || is center (know) : mostly endnode neighbors
    if(edsquare < BSplineGridNodeMap[x][y][z]->ed2c.EDSquare) {
      BSplineGridNodeMap[x][y][z]->ed2c.EDSquare = edsquare;
    }
  }else{
    BSplineGridNodeMap[x][y][z]->ed2c.EDSquare = edsquare;
    BSplineGridNodeMap[x][y][z]->ed2c.has_center = 2;
  }

  return std::max(0.0,sqrt(BSplineGridNodeMap[x][y][z]->ed2c.EDSquare));
}

double bnuk_search::getSteplen_opt(int x, int y, int z) {
  int startds;
  if (BSplineGridNodeMap[x][y][z]->ed2c.has_center == 1)// 1 visited has_center
  {
    Vector3i center = BSplineGridNodeMap[x][y][z]->ed2c.center;
    double R = sqrt(BSplineGridNodeMap[center(0)][center(1)][center(2)]->ed2c.EDSquare);
    double d = sqrt(distanceSquare_2c(x - center(0), y - center(1), z - center(2)));
    startds = (int) std::max(floor((R - d) * (R - d)) - 1.0, 0.0);
  } else if (BSplineGridNodeMap[x][y][z]->ed2c.has_center == 2)//2 visited || is center (know)
    return sqrt(BSplineGridNodeMap[x][y][z]->ed2c.EDSquare);
  else//BSplineGridNodeMap[x][y][z]->ed2c.has_center == 0  not visited  ||a new node
    startds = 0;
  bool collision_detected = false;

  int ptx, pty, ptz;
  for (int ds = startds; ds < distance_queue.size(); ds++) {
    if (!distance_queue[ds].empty()) {
      for (auto &it:distance_queue[ds]) {
        ptx = x + (it)(0);
        pty = y + (it)(1);
        ptz = z + (it)(2);

        if (ptx < 0 || ptx >= x_size || pty < 0 || pty >= y_size || ptz < 0 || ptz >= z_size) {//detected bounding
          collision_detected = true;
          BSplineGridNodeMap[x][y][z]->ed2c.EDSquare = ds;
          BSplineGridNodeMap[x][y][z]->ed2c.has_center = 2;
          BSplineGridNodeMap[x][y][z]->ed2c.center = Vector3i(x, y, z);
          break;
        }

        if (checkPtCollision(ptx, pty, ptz)) {
          collision_detected = true;
          BSplineGridNodeMap[x][y][z]->ed2c.EDSquare = ds;
          BSplineGridNodeMap[x][y][z]->ed2c.has_center = 2;
          BSplineGridNodeMap[x][y][z]->ed2c.center = Vector3i(x, y, z);
          break;
        } else //opt
        {
          if (BSplineGridNodeMap[ptx][pty][ptz]->ed2c.has_center != 2) {
            BSplineGridNodeMap[ptx][pty][ptz]->ed2c.has_center = 1;
            BSplineGridNodeMap[ptx][pty][ptz]->ed2c.center = Vector3i(x, y, z);
          }

        }
      }
    }
    if (collision_detected)
      break;
  }
  if (!collision_detected) {
    BSplineGridNodeMap[x][y][z]->ed2c.EDSquare = distance_queue.size() - 1;
    BSplineGridNodeMap[x][y][z]->ed2c.center = Vector3i(x, y, z);
  }
//    if(z>=4&&z<=20)
//    {
//        double velmy=sqrt(BSplineGridNodeMap[x][y][z]->ed2c.EDSquare);
//        double veledt=BSplineGridNodeMap[x][y][z]->steplen;
//        if(velmy!=veledt)
//        {
//            cout<<"index:"<<x<<" "<<y<<" "<<z<<endl;
//            cout<<"vel my:"<<velmy<<endl;
//            cout<<"vel EDT:"<<veledt<<endl<<endl;
//        }
//    }
  return sqrt(BSplineGridNodeMap[x][y][z]->ed2c.EDSquare);
}

double bnuk_search::getSteplen_edt(int x, int y, int z) {
  return BSplineGridNodeMap[x][y][z]->steplen;
}

inline bool bnuk_search::checkPtCollision(int x, int y, int z) {
  if (x < 0 || x >= x_size || y < 0 || y >= y_size || z < 0 || z >= z_size) {
    cerr << "checkPtCollision: index not in bound" << endl;
    return false;
  }

  if (BSplineGridNodeMap[x][y][z]->occupancy > 0.5)
    return true;
  else
    return false;
}


inline void bnuk_search::setMaxVel(double vel) {
  if (!gridnode_map_init) {
    cerr << "setMaxVel(): no resolution! Init gridnodemap first!" << endl;
  }
  max_vel = vel;
  max_vel_set = true;
  int maxdis = round(max_vel * inv_resolution);
  if (!distance_queue_init) {
    distance_queue_init = true;
    initDistanceQueue(maxdis);
  }
}

Vector4f bnuk_search::getDynamicsConstraint() {
  Vector4f dynamic_constraint;
  dynamic_constraint(0) = dynamics_max_vel;
  dynamic_constraint(1) = dynamics_min_vel;
  dynamic_constraint(2) = dynamics_max_acc;
  dynamic_constraint(3) = dynamics_min_acc;
  return dynamic_constraint;
}

void bnuk_search::setDynamicsConstraint(Vector4f dynamic_constraint) {
  dynamics_max_vel = dynamic_constraint(0);
  dynamics_min_vel = dynamic_constraint(1);
  dynamics_max_acc = dynamic_constraint(2);
  dynamics_min_acc = dynamic_constraint(3);
}
/*******other*******/

uint64_t bnuk_search::C_n_k(uint64_t n, uint64_t k) {
  if (k > n) {
    return 0;
  }
  uint64_t r = 1;
  for (uint64_t d = 1; d <= k; ++d) {
    r *= n--;
    r /= d;
  }
  return r;
}


void bnuk_search::initDistanceQueue(int maxDis) {
  int maxDisSquare = maxDis * maxDis;
  distance_queue.resize(maxDisSquare + 1);
  int dx, dy, dz;
  int dis;
  for (dx = -maxDis; dx <= maxDis; dx++) {
    for (dy = -maxDis; dy <= maxDis; dy++) {
      for (dz = -maxDis; dz <= maxDis; dz++) {
        dis = distanceSquare_2c(dx, dy, dz);
        if (dis <= maxDisSquare)
          distance_queue[dis].push_back(Vector3i(dx, dy, dz));
      }
    }
  }
}


#endif //PROJECT_BNUK_SEARCH_IMPL