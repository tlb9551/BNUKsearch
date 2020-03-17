
//
// Created by tlb on 19-7-27.
//

#ifndef BNUK_OPT_B_SPLINE_OPTIMIZER_H
#define BNUK_OPT_B_SPLINE_OPTIMIZER_H
#include "mosek.h"
#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
class matrix_helper {
  friend class b_spline_optimizer;
 private:
  double dt = 1.0;
  double pow_inv_dt_[2 * 6];

  static Eigen::Matrix<double, 6, 6> bspline_basis_matrix;
  static Eigen::Matrix<double, 6, 6> base_coefficients;
  Eigen::Matrix<double, 6, 6> diff_matrix;

  static uint64_t C_n_k(uint64_t n, uint64_t k) {
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
  static Eigen::Matrix<double, 6, 6> computeBlendingMatrix() {
    Eigen::Matrix<double, 6, 6> m;
    m.setZero();

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double sum = 0;

        for (int s = j; s < 6; ++s) {
          sum += std::pow(-1.0, s - j) * C_n_k(6, s - j)
              * std::pow(6 - s - 1.0, 6 - 1.0 - i);
        }
        m(i, j) = C_n_k(6 - 1, 6 - 1 - i) * sum;
      }
    }

    uint64_t factorial = 1.0;
    for (int i = 2; i < 6; ++i) {
      factorial *= i;
    }

    return (m / factorial).template cast<double>();
  }
  static Eigen::Matrix<double, 6, 6> computeBaseCoefficients() {
    Eigen::Matrix<double, 6, 6> base_coefficients_;
    base_coefficients_.setZero();
    base_coefficients_.row(0).setOnes();

    const int DEG = 5;
    int order = DEG;
    for (int n = 1; n < DEG + 1; n++) {
      for (int i = DEG - order; i < DEG + 1; i++) {
        base_coefficients_(n, i) = (order - DEG + i) * base_coefficients_(n - 1, i);
      }
      order--;
    }
    return base_coefficients_.template cast<double>();
  }
  /*  Q matrix without 1/dt^n **/
  static std::array<Eigen::Matrix<double, 6, 6>, 4> computeQuadraticCoefficients() {
    std::array<Eigen::Matrix<double, 6, 6>, 4> res;

    const int DEG = 5;

    for (int derivative = 1; derivative < 5; derivative++) {
      res[derivative - 1].setZero();

      for (int col = 0; col < 6 - derivative; col++) {
        for (int row = 0; row < 6 - derivative; row++) {
          double exp = (DEG - derivative) * 2 + 1 - row - col;

          res[derivative - 1](DEG - row, DEG - col) =
              base_coefficients(derivative, 6 - 1 - row) *
                  base_coefficients(derivative, 6 - 1 - col) * 2 / exp;
        }
      }
    }

    return res;
  }
  /*  M^T*Q*M matrix with 1/dt^n **/
  std::array<Eigen::Matrix<double, 6, 6>, 4> computeQuadraticCostJacobian() {
    std::array<Eigen::Matrix<double, 6, 6>, 4> m_array;

    std::array<Eigen::Matrix<double, 6, 6>, 4>
        quadratic_coefficients = computeQuadraticCoefficients();

    for (int i = 0; i < 4; ++i) {
      Eigen::Matrix<double, 6, 6> m;
      m = pow_inv_dt_[2 * (i + 1) - 1] * bspline_basis_matrix.transpose() *
          quadratic_coefficients[i] * bspline_basis_matrix;

      m_array[i] = m.template cast<double>();
    }
    return m_array;
  }
  static Eigen::Matrix<double, 6, 6> computeBasisMatrix_Bezier(int diff_degree) {
    if (diff_degree > 5) {
      std::cerr << "computeBasisMatrix_Bezier: diff_degree  = " << diff_degree
                << std::endl;
      return Eigen::Matrix<double, 6, 6>::Zero();
    }

    Eigen::Matrix<double, 6, 6> m;
    m.setZero();

    for (int i = 0; i < 5 + 1 - diff_degree; ++i) {
      for (int j = 0; j < 5 + 1 - diff_degree; ++j) {
        if (i < j)
          m(i, j) = 0;
        else {
          m(i, j) = std::pow(-1.0, i - j) * C_n_k(5 - diff_degree, j) * C_n_k(5 - diff_degree - j, i - j);
        }
      }
    }

    if (diff_degree > 0) {
      for (int i = 0; i < diff_degree; i++)
        m(5 - i, 5 - i) = 1.0;
    }

    return (m).template cast<double>();
  }

  explicit matrix_helper(double _dt) {
    dt = _dt;

    /**init pow_inv_dt_  **/
    double current_pow_inv_dt = 1.0;
    for (int i = 0; i < 2 * 6; ++i) {
      pow_inv_dt_[i] = current_pow_inv_dt;
      current_pow_inv_dt /= dt;
    }
    /**init diff_matrix  **/
    diff_matrix = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 5; i++)
      diff_matrix(i, i + 1) = i + 1;
    diff_matrix = diff_matrix / dt;

  };

  Eigen::Matrix<double, 6, 6> inputDiff0() {
    Eigen::Matrix<double, 6, 6> res = computeBasisMatrix_Bezier(0).inverse() * bspline_basis_matrix;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        if(res(i,j) < 1e-9 && res(i,j) > -1e-9)
          res(i,j)=0;
      }
    }
    return res;
  }
  Eigen::Matrix<double, 6, 6> inputDiff1() {
    Eigen::Matrix<double, 6, 6> res = computeBasisMatrix_Bezier(1).inverse() * diff_matrix * bspline_basis_matrix;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        if(res(i,j) < 1e-9 && res(i,j) > -1e-9)
          res(i,j)=0;
      }
    }
    return res;
  }
  Eigen::Matrix<double, 6, 6> inputDiff2() {
    Eigen::Matrix<double, 6, 6> res = computeBasisMatrix_Bezier(2).inverse() * diff_matrix * diff_matrix * bspline_basis_matrix;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        if(res(i,j) < 1e-9 && res(i,j) > -1e-9)
          res(i,j)=0;
      }
    }
    return res;
  }
  Eigen::Matrix<double, 5, 5> inputEndMatrix() {
    Eigen::Matrix<double, 5 + 1, 5 + 1> endm;
    for (int i = 0; i < 5 + 1; i++) {
      for (int j = 0; j < 5 + 1; j++) {
        if (j >= i) {
          int factorial = 1;
          for (int k = 0; k < i; k++)
            factorial *= (j - k);
          endm(i, j) = factorial * 1.0 / pow_inv_dt_[j];
        } else endm(i, j) = 0.0;
      }
    }
    endm = endm * bspline_basis_matrix;
    return endm.block(0, 1, 5, 5);
  }

};
class b_spline_optimizer {
 public:
  int DEG = 5;
  double dt = 1.0;
  /***mode***/
  bool enableQCon = false;
  /***param***/
  double max_vel = 3.0;
  double max_acc = 3.0;
  Eigen::Vector3d map_ld_corner;
  Eigen::Vector3d map_ru_corner;

  Eigen::Vector3d end_position;
  Eigen::Vector3d end_velocity;

  /***var***/
  double *xx;
  /***Q_obj***/
  MSKint32t *qsubi;
  MSKint32t *qsubj;
  double *qval;
  int num_qnz;
  int idx_q = 0;
  double *weight;
  std::array<Eigen::Matrix<double, 6, 6>, 4> MQM_array;
  /***C***/
  MSKint32t *csubi;
  double *cval;
  /***Cf***/
  MSKrealt cf = 0.0;
  /***A***/
  MSKint32t *asubi;
  MSKint32t *asubj;
  double *aval;
  int idx_a = 0;
  int arows = 0;
  int num_anz;
  int num_cons;
  /***Q_con***/
  MSKint32t *qcsubi;
  MSKint32t *qcsubj;
  double *qcval;
  int num_q_con = 0;
  /* Bounds on constraints. */
  MSKboundkeye *bkc;
  double *blc;
  double *buc;

  MSKboundkeye *bkx;
  double *blx;
  double *bux;

  /***Var***/
  int numvar_dim;
  int numvar;

  /***Input Spline***/
  std::vector<Eigen::Vector3d> control_pts;

  /***Matrix***/
  Eigen::Matrix<double, 6, 6> MQM;
  Eigen::Matrix<double, 6, 6> bspline2bezier_diff0;
  Eigen::Matrix<double, 6, 6> bspline2bezier_diff1;
  Eigen::Matrix<double, 6, 6> bspline2bezier_diff2;
  Eigen::Matrix<double, 5, 5> end_matrix;

  /***func***/
  bool solve();

  void setControlPts(const std::vector<Eigen::Vector3d> &pts);
  void setWeight(const std::vector<double> &w);
  void setMQM(int idx = 2);
  void setQobjVal();
  void setCobjVal_SafeLimit();
  void setAVal_DynamicCons();
  void setAVal_EndCons();
  void setQconVal(std::vector<Eigen::Vector3d> centers,std::vector<double> radius);
  void usingQcon();
  void setEndPos(Eigen::Vector3d pos) { end_position = pos; };
  void setEndVel(Eigen::Vector3d vel) { end_velocity = vel; };
  void setMapRange(Eigen::Vector3d ld,Eigen::Vector3d ru){map_ld_corner = ld;map_ru_corner = ru;};
  void setDynamicCons(double _max_vel,double _max_acc){max_vel = _max_vel;max_acc = _max_acc;};
  void printProblem();
  b_spline_optimizer() : helper(dt) {
    setMQM();
    MQM_array = helper.computeQuadraticCostJacobian();
    end_matrix = helper.inputEndMatrix();
    bspline2bezier_diff0 = helper.inputDiff0();
    bspline2bezier_diff1 = helper.inputDiff1();
    bspline2bezier_diff2 = helper.inputDiff2();
  };
  b_spline_optimizer(double dt_) : helper(dt) {
    dt = dt_;
    setMQM();
    MQM_array = helper.computeQuadraticCostJacobian();
    end_matrix = helper.inputEndMatrix();
    bspline2bezier_diff0 = helper.inputDiff0();
    bspline2bezier_diff1 = helper.inputDiff1();
    bspline2bezier_diff2 = helper.inputDiff2();

  };
  ~b_spline_optimizer() {
    delete[] xx;
    delete[] qsubi;
    delete[] qsubj;
    delete[] qval;
    delete[] weight;
    delete[] csubi;
    delete[] cval;
    delete[] asubi;
    delete[] asubj;
    delete[] aval;
    delete[] bkc;
    delete[] blc;
    delete[] buc;
    delete[] bkx;
    delete[] blx;
    delete[] bux;
    delete[] qcsubi;
    delete[] qcsubj;
    delete[] qcval;
  }
 private:
  matrix_helper helper;

};

Eigen::Matrix<double, 6, 6> matrix_helper::bspline_basis_matrix = computeBlendingMatrix();
Eigen::Matrix<double, 6, 6> matrix_helper::base_coefficients = computeBaseCoefficients();

#endif //BNUK_OPT_B_SPLINE_OPTIMIZER_H
