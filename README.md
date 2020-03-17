# BNUKsearch
## 1.What's this?

BNUKsearch is a local trajectory planner for quadrotor trajectory planning. In this ROS package we show an example of trajectory planning and navigation using the BNUK algorithm under the RHC framework. This is only a research code.

For details we refer readers to our paper.

* **Real-time Trajectory Generation for Quadrotors using B-spline based Non-uniform Kinodynamic Search,** Lvbang Tang, Hesheng Wang, Peng Li, Yong Wang, 2019 IEEE International Conference on Robotics and Biomimetics (ROBIO)

If you use this planning algorithm for your academic research, please cite our paper.
```
@inproceedings{tang2019real,
  title={Real-time Trajectory Generation for Quadrotors using B-spline based Non-uniform Kinodynamic Search},
  author={Tang, Lvbang and Wang, Hesheng and Li, Peng and Wang, Yong},
  booktitle={2019 IEEE International Conference on Robotics and Biomimetics (ROBIO)},
  pages={1133--1138},
  year={2019},
  organization={IEEE}
}
```
## 2.How to install?

Our testing environment: **Ubuntu** 16.04, **ROS** Kinetic.

- Step 1: You should install the [armadillo](http://arma.sourceforge.net/), which is a c++ linear algebra library. Then clone and compile the modified [plan_utils](https://github.com/tlb9551/plan_utils), which contains several ROS-package used for running the simulation.
```
  sudo apt-get install libarmadillo-dev
  cd ~/catkin_ws/src
  git clone https://github.com/tlb9551/plan_utils.git
  cd ..
  catkin_make
  source ~/catkin_ws/devel/setup.bash
```

- Step 2: Clone the repository to your catkin workspace and catkin_make. For example:
```
  cd ~/catkin_ws/src
  git clone https://github.com/tlb9551/BNUKsearch.git
  cd ..
  catkin_make
  source ~/catkin_ws/devel/setup.bash
```

- Step 3: Install **Mosek** license: We use **Mosek** for solving quadratically constrainted quadratic programming (QCQP). An [academic license](https://www.mosek.com/products/academic-licenses/) is required. Then put the license file (**.lic) into ```$HOME/mosek/mosek.lic```.

#### common issues:

- 

## 3.How to use?
After completing the above installation and compilation, you can run the navigation simulation example with the follows
```
  roslaunch bnuk_opt bnuk_replan.launch
```
In rviz, you can see a plugin 'Goal3DTool' in toolbar. It is used to publish a goal position. To use it, click the 'Goal3DTool' (or press keyboard 'g'), then press on left mouse button on a position in rviz, click right mouse button to start to drag it slide up or down for a targeting height (don't loose left button at this time). Finally you loose left mouse button and a goal position will be published.

## 4.Acknowledgements
This simulation framework was inspired by [Fei Gao](https://ustfei.com/)'s [Btraj](https://github.com/HKUST-Aerial-Robotics/Btraj), the basic framework was created by him and is open source, thanks for his excellent work and selflessness! I believe that open source can make research in this field more active!

We also use [mosek](https://www.mosek.com/) for solving quadratically constrainted quadratic programming(QCQP), [sdf_tools](https://github.com/UM-ARM-Lab/sdf_tools) for building euclidean distance field.

## 5.Licence
The source code is released under [GPLv3](http://www.gnu.org/licenses/) license.

## 6.Notes
This code is only used for research work. It is not fully optimized. If you find any problems or have any ideas for improvement, you can raise a issue.
