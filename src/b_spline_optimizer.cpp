//
// Created by tlb on 19-7-27.
//

#include "b_spline_optimizer.h"
using namespace std;
static void MSKAPI printstr(void *handle, MSKCONST char str[]) {
//  printf("%s", str);
}
bool b_spline_optimizer::solve() {
  MSKrescodee r;
  MSKenv_t env;
  MSKtask_t task;
  bool res = false;
/* Create the mosek environment. */
  r = MSK_makeenv(&env, NULL);

  if (r == MSK_RES_OK) {
    /* Create the optimization task. */
    r = MSK_maketask(env, num_cons, numvar, &task);

    if (r == MSK_RES_OK) {
      r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);


      /* Append 'NUMCON' empty constraints.
       The constraints will initially have no bounds. */
      if (r == MSK_RES_OK)
        r = MSK_appendcons(task, num_cons);


      /* Append 'NUMVAR' variables.
       The variables will initially be fixed at zero (x=0). */
      if (r == MSK_RES_OK)
        r = MSK_appendvars(task, numvar);

      /* Optionally add a constant term to the objective. */
      if (r == MSK_RES_OK)
        r = MSK_putcfix(task, cf);

      /* Set the linear term c_j in the objective.*/
      if (r == MSK_RES_OK)
        r = MSK_putclist(task, numvar, csubi, cval);

      for (int j = 0; j < numvar && r == MSK_RES_OK; ++j) {
        /* Set the bounds on variable j.
         blx[j] <= x_j <= bux[j] */
        if (r == MSK_RES_OK)
          r = MSK_putvarbound(task,
                              j,           /* Index of variable.*/
                              bkx[j],      /* Bound key.*/
                              blx[j],      /* Numerical value of lower bound.*/
                              bux[j]);     /* Numerical value of upper bound.*/
      }
      /* Input A */
      if (r == MSK_RES_OK)
        r = MSK_putaijlist(task,
                           num_anz,
                           asubi,
                           asubj,
                           aval);

      /* Set the bounds on constraints.
         for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
      for (int i = 0; i < num_cons && r == MSK_RES_OK; ++i)
        r = MSK_putconbound(task,
                            i,           /* Index of constraint.*/
                            bkc[i],      /* Bound key.*/
                            blc[i],      /* Numerical value of lower bound.*/
                            buc[i]);     /* Numerical value of upper bound.*/

      if (enableQCon) {
        /* Input the Q  for the cons. */

        for (int k = 0; k < num_q_con; ++k) {
          if (r == MSK_RES_OK) {
            r = MSK_putqconk(task,
                             k + num_cons - num_q_con,
                             21 * 3,
                             &qcsubi[21 * 3 * k],
                             &qcsubj[21 * 3 * k],
                             &qcval[21 * 3 * k]);
          }

        }

      }

      /* Input the Q for the objective. */
      if (r == MSK_RES_OK)
        r = MSK_putqobj(task, num_qnz, qsubi, qsubj, qval);

      if (r == MSK_RES_OK)
        r = MSK_putparam(task, "MSK_IPAR_INTPNT_MAX_ITERATIONS", "400");
      if (r == MSK_RES_OK) {

      } else {
        std::cerr << "error in param setting" << std::endl;
      }

      if (r == MSK_RES_OK) {
        MSKrescodee trmcode;

        /* Run optimizer */
        r = MSK_optimizetrm(task, &trmcode);

        /* Print a summary containing information
           about the solution for debugging purposes*/
        MSK_solutionsummary(task, MSK_STREAM_MSG);

        if (r == MSK_RES_OK) {
          MSKsolstae solsta;

          MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

          switch (solsta) {
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
              MSK_getxx(task,
                        MSK_SOL_ITR,    /* Request the interior solution. */
                        xx);

              printf("Optimal primal solution\n");
              res = true;
//              for (int j = 0; j < numvar; ++j)
//                printf("x[%d]: %e\n", j, xx[j]);

              break;

            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:printf("Primal or dual infeasibility certificate found.\n");
              break;

            case MSK_SOL_STA_UNKNOWN:printf("The status of the solution could not be determined.\n");
              break;

            default:printf("Other solution status.");
              break;
          }
        } else {
          printf("Error while optimizing.\n");
        }
      }

      if (r != MSK_RES_OK) {
        /* In case of an error print error code and description. */
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];

        printf("An error occurred while optimizing.\n");
        MSK_getcodedesc(r,
                        symname,
                        desc);
        printf("Error %s - '%s'\n", symname, desc);
      }
    }
    MSK_deletetask(&task);
  }
  MSK_deleteenv(&env);
  return (res);
}

void b_spline_optimizer::setControlPts(const std::vector<Eigen::Vector3d> &pts) {
  control_pts = pts;
  numvar_dim = control_pts.size();
  numvar = numvar_dim * 3;
  xx = new double[numvar];

  num_qnz = 3 * ((2 * numvar_dim - DEG) * (DEG + 1) / 2);
  qsubi = new MSKint32t[num_qnz];
  qsubj = new MSKint32t[num_qnz];
  qval = new double[num_qnz];
  weight = new double[numvar_dim];
  idx_q = 0;
  csubi = new MSKint32t[numvar];
  cval = new double[numvar];

  int NUMANZ_vel = 20 * (numvar_dim - DEG) * 3;
  int NUMANZ_acc = 16 * (numvar_dim - DEG) * 3;
  int NUMANZ_end = DEG * 2 * 3;
  int NUMCONS_vel = 5 * (numvar_dim - DEG) * 3;
  int NUMCONS_acc = 4 * (numvar_dim - DEG) * 3;
  int NUMCONS_end = 2 * 3;
  int NUMCONS_start = DEG * 3;
  int NUMCONS_Qcon = (numvar_dim - DEG) * (DEG + 1);
  int NUMANZ_Qcon = NUMCONS_Qcon * (DEG + 1) * 3;
  if (enableQCon) {
    num_anz = NUMANZ_acc + NUMANZ_vel + NUMANZ_end + NUMANZ_Qcon;
    num_cons = NUMCONS_vel + NUMCONS_acc + NUMCONS_end + NUMCONS_Qcon;
    num_q_con = NUMCONS_Qcon;
  } else {
    num_anz = NUMANZ_acc + NUMANZ_vel + NUMANZ_end;
    num_cons = NUMCONS_vel + NUMCONS_acc + NUMCONS_end;
  }
  qcsubi = new MSKint32t[21 * 3 * num_q_con];//lower triangular part of Qc is inputted, so 21 num vals per dim
  qcsubj = new MSKint32t[21 * 3 * num_q_con];
  qcval = new double[21 * 3 * num_q_con];

  asubi = new MSKint32t[num_anz];
  asubj = new MSKint32t[num_anz];
  aval = new double[num_anz];
  bkc = new MSKboundkeye[num_cons];
  blc = new double[num_cons];
  buc = new double[num_cons];

  bkx = new MSKboundkeye[numvar];
  blx = new double[numvar];
  bux = new double[numvar];

  for (int i = 0; i < DEG; ++i) {
    for (int dim = 0; dim < 3; ++dim) {
      bkx[i + dim * numvar_dim] = MSK_BK_FX;
      blx[i + dim * numvar_dim] = control_pts[i](dim);
      bux[i + dim * numvar_dim] = control_pts[i](dim);
    }
  }
  for (int i = DEG; i < numvar_dim; ++i) {
    for (int dim = 0; dim < 3; ++dim) {
      bkx[i + dim * numvar_dim] = MSK_BK_RA;
      blx[i + dim * numvar_dim] = map_ld_corner(dim);
      bux[i + dim * numvar_dim] = map_ru_corner(dim);
    }
  }

  idx_a = 0;
  arows = 0;
}
void b_spline_optimizer::setWeight(const std::vector<double> &w) {
  assert(weight != nullptr);
  assert(w.size() == numvar_dim);
  for (int i = 0; i < w.size(); i++) {
    weight[i] = w[i];
  }
}

void b_spline_optimizer::setMQM(int idx) {
  MQM = MQM_array[idx];
}

void b_spline_optimizer::setQobjVal() {
  double sumtrace[DEG + 1];
  for (int i = 0; i < DEG + 1; i++) {
    sumtrace[i] = MQM.diagonal(i).sum();
  }

  int row_offset = 0;
  int col_offset = 0;
  for (int dim = 0; dim < 3; ++dim) {
    row_offset = dim * numvar_dim;
    col_offset = dim * numvar_dim;
    for (int i = 1; i < DEG + 1; i++) {
      for (int dx = 0; dx < numvar_dim - i; dx++) {
        if (dx >= DEG - i && dx <= numvar_dim - (DEG - i) - i - 1) {
//          qsubi[idx_q] = dx + row_offset;
//          qsubj[idx_q] = dx + i + col_offset;
//          qval[idx_q] = 2*sumtrace[i];
//          idx_q++;

          qsubi[idx_q] = dx + i + row_offset;
          qsubj[idx_q] = dx + col_offset;
          qval[idx_q] = 2 * sumtrace[i];
          idx_q++;
        } else {
          auto diag = MQM.diagonal(i);
          int abs_dx = std::min(dx, abs(numvar_dim - i - dx - 1));
          double diag_part_sum = 0;
          for (int k = 0; k <= abs_dx; k++) {
            diag_part_sum += diag(k);
          }
//          qsubi[idx_q] = dx + row_offset;
//          qsubj[idx_q] = dx + i + col_offset;
//          qval[idx_q] = 2*diag_part_sum;
//          idx_q++;

          qsubi[idx_q] = dx + i + row_offset;
          qsubj[idx_q] = dx + col_offset;
          qval[idx_q] = 2 * diag_part_sum;
          idx_q++;
        }
      }
    }
    {//diag
      int i = 0;
      for (int dx = 0; dx < numvar_dim - i; dx++) {
        if (dx >= DEG - i && dx <= numvar_dim - (DEG - i) - i - 1) {
          qsubi[idx_q] = dx + row_offset;
          qsubj[idx_q] = dx + i + col_offset;
          qval[idx_q] = 2 * (sumtrace[i] + weight[i]);
          idx_q++;
        } else {
          auto diag = MQM.diagonal(i);
          int abs_dx = std::min(dx, abs(numvar_dim - i - dx - 1));
          double diag_part_sum = 0;
          for (int k = 0; k <= abs_dx; k++) {
            diag_part_sum += diag(k);
          }
          qsubi[idx_q] = dx + row_offset;
          qsubj[idx_q] = dx + i + col_offset;
          qval[idx_q] = 2 * (diag_part_sum + weight[i]);
          idx_q++;
        }
      }
    }
  }

}
void b_spline_optimizer::setCobjVal_SafeLimit() {
  assert(weight != nullptr);
  assert(numvar_dim != 0);
  assert(!control_pts.empty());

  int idx = 0;
  for (int dim = 0; dim < 3; ++dim) {
    for (int i = 0; i < numvar_dim; ++i) {
      csubi[idx] = idx;
      cval[idx] = -weight[i] * 2 * control_pts[i](dim);
      cf += weight[i] * control_pts[i](dim) * control_pts[i](dim);
      idx++;
    }
  }

}

void b_spline_optimizer::setAVal_DynamicCons() {

  double EPS = 0.00001;
  Eigen::Matrix<double, 5, 6> b2b_vel = bspline2bezier_diff1.block(0, 0, 5, 6);
  Eigen::Matrix<double, 4, 6> b2b_acc = bspline2bezier_diff2.block(0, 0, 4, 6);

  MSKint32t b2b_vel_subi[20];
  MSKint32t b2b_vel_subj[20];
  double b2b_vel_val[20];
  int b2b_vel_nonzero = 0;
  for (int row = 0; row < 5; ++row) {
    for (int col = 0; col < 6; ++col) {
      if (b2b_vel(row, col) > EPS || b2b_vel(row, col) < -EPS) {
        b2b_vel_subi[b2b_vel_nonzero] = row;
        b2b_vel_subj[b2b_vel_nonzero] = col;
        b2b_vel_val[b2b_vel_nonzero] = b2b_vel(row, col);
        b2b_vel_nonzero++;
      }
    }
  }

  MSKint32t b2b_acc_subi[16];
  MSKint32t b2b_acc_subj[16];
  double b2b_acc_val[16];
  int b2b_acc_nonzero = 0;
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 6; ++col) {
      if (b2b_acc(row, col) > EPS || b2b_acc(row, col) < -EPS) {
        b2b_acc_subi[b2b_acc_nonzero] = row;
        b2b_acc_subj[b2b_acc_nonzero] = col;
        b2b_acc_val[b2b_acc_nonzero] = b2b_acc(row, col);
        b2b_acc_nonzero++;
      }
    }
  }

  int row_offset = arows;
  int col_offset = 0;
  for (int dim = 0; dim < 3; ++dim) {

    //vel part
    row_offset = arows + dim * (5 * (numvar_dim - DEG) + 4 * (numvar_dim - DEG));
    col_offset = dim * numvar_dim;
    for (int block_num = 0; block_num < numvar_dim - DEG; ++block_num) {
      for (int i = 0; i < b2b_vel_nonzero; i++) {
        asubi[idx_a] = b2b_vel_subi[i] + row_offset;
        asubj[idx_a] = b2b_vel_subj[i] + col_offset;
        aval[idx_a] = b2b_vel_val[i];
        idx_a++;
      }
      for (int i = 0; i < 5; i++) {
        bkc[i + row_offset] = MSK_BK_RA;
        blc[i + row_offset] = -max_vel;
        buc[i + row_offset] = +max_vel;
      }
      row_offset += 5;
      col_offset += 1;
    }

    //acc part
    col_offset = dim * numvar_dim;
    for (int block_num = 0; block_num < numvar_dim - DEG; ++block_num) {
      for (int i = 0; i < b2b_acc_nonzero; i++) {
        asubi[idx_a] = b2b_acc_subi[i] + row_offset;
        asubj[idx_a] = b2b_acc_subj[i] + col_offset;
        aval[idx_a] = b2b_acc_val[i];
        idx_a++;
      }
      for (int i = 0; i < 4; i++) {
        bkc[i + row_offset] = MSK_BK_RA;
        blc[i + row_offset] = -max_acc;
        buc[i + row_offset] = +max_acc;
      }
      row_offset += 4;
      col_offset += 1;
    }

  }
  arows = row_offset;

}
void b_spline_optimizer::setAVal_EndCons() {
  int row_offset = arows;
  for (int dim = 0; dim < 3; ++dim) {
    for (int col = 0; col < 5; col++) {//end position
      asubi[idx_a] = row_offset;
      asubj[idx_a] = (dim + 1) * numvar_dim - 5 + col;
      aval[idx_a] = end_matrix(0, col);
      idx_a++;
    }
    bkc[row_offset] = MSK_BK_FX;
    blc[row_offset] = end_position(dim);
    buc[row_offset] = end_position(dim);
    row_offset++;

    for (int col = 0; col < 5; col++) {//end velocity
      asubi[idx_a] = row_offset;
      asubj[idx_a] = (dim + 1) * numvar_dim - 5 + col;
      aval[idx_a] = end_matrix(1, col);
      idx_a++;
    }
    bkc[row_offset] = MSK_BK_FX;
    blc[row_offset] = end_velocity(dim);
    buc[row_offset] = end_velocity(dim);
    row_offset++;
  }
  arows = row_offset;
}
void b_spline_optimizer::setQconVal(std::vector<Eigen::Vector3d> centers, std::vector<double> radius) {
  if (enableQCon) {
    if (radius.size() != control_pts.size() - DEG) {
      std::cerr << "radius or centers's size not correct: " << centers.size() << " and " << radius.size() << std::endl;
      return;
    }

    std::vector<std::vector<MSKint32t>> Q_subi(6);
    std::vector<std::vector<MSKint32t>> Q_subj(6);
    std::vector<std::vector<double>> Q_val(6);

    for (int i = 0; i < DEG + 1; i++) {
      Eigen::Matrix<double, 6, 6> Q;
      Q = bspline2bezier_diff0.block(i, 0, 1, DEG + 1).transpose()
          * bspline2bezier_diff0.block(i, 0, 1, DEG + 1);
      for (int row = 0; row < DEG + 1; row++) {
        for (int col = 0; col < DEG + 1; ++col) {
          if (row >= col) {
            Q_subi[i].push_back(row);
            Q_subj[i].push_back(col);
            Q_val[i].push_back(Q(row, col));
          }
        }
      }
    }

    int idx_Q = 0;
    int idx_A = idx_a;
    int row_offset = arows;
    for (int i = 0; i < numvar_dim - DEG; ++i) {
      double r = radius[i];
      for (int j = 0; j < DEG + 1; ++j) {
        for (int dim = 0; dim < 3; ++dim) {
          for (int idx = 0; idx < Q_subi[j].size(); idx++) {
            qcsubi[idx_Q] = Q_subi[j][idx] + dim * numvar_dim + i;
            qcsubj[idx_Q] = Q_subj[j][idx] + dim * numvar_dim + i;
            qcval[idx_Q] = 2 * Q_val[j][idx];
            idx_Q++;
          }
          for (int local_cps = 0; local_cps < DEG + 1; ++local_cps) {
            asubi[idx_A] = row_offset;
            asubj[idx_A] = i + local_cps + dim * numvar_dim;
            aval[idx_A] = -2.0 * centers[i](dim) * bspline2bezier_diff0(j, local_cps);
            idx_A++;
          }
        }
        bkc[row_offset] = MSK_BK_UP;
        blc[row_offset] = -MSK_INFINITY;
        buc[row_offset] = r * r - centers[i].squaredNorm();
        row_offset++;
      }
    }

    idx_a = idx_A;
    arows = row_offset;
    assert(arows == num_cons);
    assert(idx_a == num_anz);
  }
}
void b_spline_optimizer::usingQcon() {
  enableQCon = true;
}
void b_spline_optimizer::printProblem() {

  cout << "Num of constraints:" << num_cons << endl;
  cout << "Num of constraints:" << numvar << endl;
  cout << "Constant term of the objective:" << cf << endl;
  cout << "Linear term c_j of the objective:" << endl;
  for (int i = 0; i < numvar; i++) {
    cout << cval[i] << " ";
  }
  cout << endl;
  cout << "The bounds on variable:" << endl;
  for (int i = 0; i < numvar; i++) {
    switch (bkx[i]) {
      case MSK_BK_FX:cout << blx[i] << "\t=" << "\tx" << i;break;
      case MSK_BK_FR:cout << "free" << "\t=" << "\tx" << i;break;
      case MSK_BK_LO:cout << blx[i] << "\t<=" << "\tx" << i;break;
      case MSK_BK_RA:cout << blx[i] << "\t<=" << "\tx" << i << "\t<=" << bux[i];break;
      case MSK_BK_UP:cout << "\tx" << i << "\t<=" << bux[i];break;
    }
    cout<<endl;
  }
  cout<<endl;
  cout << "The A cons matrix:" << endl;
  unordered_map<int,double> hash_A;
  for(int i =0;i<num_anz;i++) {
    hash_A.insert({asubi[i]*numvar+ asubj[i], aval[i]});
  }
  for(int i = 0;i<num_cons;i++) {
    for (int j = 0; j < numvar; j++) {
      auto ptr = hash_A.find(i*numvar+j);
      if(ptr!=hash_A.end()){
        cout << ptr->second << " ";
      }else{
        cout <<0<< " ";
      }
    }
    cout << endl;
  }
  cout<<endl;
  cout<<"The bounds on constraints"<<endl;
  for (int i = 0; i < num_cons; i++) {
    switch (bkc[i]) {
      case MSK_BK_FX:cout << blc[i] <<"\t=" << "\t" << i;break;
      case MSK_BK_FR:cout << "free" << "\t=" << "\t" << i;break;
      case MSK_BK_LO:cout << blc[i] << "\t<=" << "\t" << i;break;
      case MSK_BK_RA:cout << blc[i] << "\t<=" << "\t" << i << "\t<=" << buc[i];break;
      case MSK_BK_UP:cout << "\t" << i << "\t<=" << buc[i];break;
    }
    cout<<endl;
  }
  cout<<endl;
  cout<<"The Q for the objective"<<endl;
  unordered_map<int,double> hash_Q_o;
  for(int i =0;i<num_qnz;i++) {
    hash_Q_o.insert({qsubi[i]*numvar+ qsubj[i], qval[i]});
  }
  for(int i = 0;i<numvar;i++) {
    for (int j = 0; j < numvar; j++) {
      auto ptr = hash_Q_o.find(i*numvar+j);
      if(ptr!=hash_Q_o.end()){
        cout << ptr->second << " ";
      }else{
        cout <<0<< " ";
      }
    }
    cout << endl;
  }
  cout<<endl;
//  if (enableQCon) {
//    cout << "The Q for the constraints" << endl;
//    unordered_map<int, double> hash_Q_c;
//    for (int i = 0; i < num_qnz; i++) {
//      hash_Q_o.insert({qsubi[i] * numvar + qsubj[i], qval[i]});
//    }
//    for (int i = 0; i < numvar; i++) {
//      for (int j = 0; j < numvar; j++) {
//        auto ptr = hash_Q_o.find(i * numvar + j);
//        if (ptr != hash_Q_o.end()) {
//          cout << ptr->second << " ";
//        } else {
//          cout << 0 << " ";
//        }
//      }
//      cout << endl;
//    }
//    cout << endl;
//
//  }
//    /* Input the Q  for the cons. */
//
//    for (int k = 0; k < num_q_con; ++k) {
//      if (r == MSK_RES_OK) {
//        r = MSK_putqconk(task,
//                         k + num_cons - num_q_con,
//                         21 * 3,
//                         &qcsubi[21 * 3 * k],
//                         &qcsubj[21 * 3 * k],
//                         &qcval[21 * 3 * k]);
//      }
//
//    }
//
//  }

}
