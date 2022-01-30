#include <uWS/uWS.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"
#include "helpers.h"
#include "json.hpp"

#include "spline.h"

// for convenience
using namespace Eigen;

using nlohmann::json;
using std::string;
using std::vector;
/*
这种方式定义的结构体，使用时需要先声明<结构变量>
typedef struct{
  类型 变量名;
} Student;
Student x; // 声明变量
-------------------------------------------
typedef struct Node {
  int data;
  struct Node* next;
} node;
第一步：
struct Node {
  int data;
  struct Node* next;
};
创建Node结构类型;
第二步:
typedef Node node; // 把Node数据类型命名为node
-------------------------------------------
struct Student {
  类型 变量名;
} 结构变量; // 声明变量
-------------------------------------------
struct Student {
  类型 变量名;
};
struct Student x; // 声明变量
*/
typedef struct Vehicle {
  int id;
  double x;
  double y;
  double vx;
  double vy;
  double s;
  double d;
  double speed;
} Vehicle;

typedef struct Planner {
  double target_d;
  vector<Vehicle> obstacles;
  Vehicle target_to_follow;
  int following_target_id;
  double dist_to_target;
  Eigen::MatrixXd s_trajectories;
  Eigen::VectorXd s_costs;
  Eigen::MatrixXd d_trajectories;
  Eigen::VectorXd d_costs;
  bool obstacle_following;
  bool feasible_traj_exist;
  int optimal_s_id;
  int optimal_d_id;
  double minimal_cost;
  int iters;
} Planner;

// get the lane id of ego car
int getMyLane(double d0) {
  int mylane = 1;
  if (d0 > 0 && d0 <= 4) {
    mylane = 0;
  } else if (d0 > 4 && d0 <= 8) {
    mylane = 1;
  } else {
    mylane = 2;
  }
  return mylane;
}

string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    std::istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  Eigen::VectorXd optimal_s_coeff(6);
  Eigen::VectorXd optimal_d_coeff(6);
  double s_cost = 999;
  double d_cost = 999;
  int step = 0;

  int minus_max_s = 0;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
               &map_waypoints_dx,&map_waypoints_dy,
               &optimal_s_coeff, &step, &max_s, &minus_max_s]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
               uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          /**
           * TODO: define a path made up of (x,y) points that the car will visit
           *   sequentially every .02 seconds
           */
          double pos_x;
          double pos_y;
          double pos_yaw;

          double MPH2mps = 1.0 / 2.23694;
          double max_speed = 43 * MPH2mps;

          int prev_path_size = previous_path_x.size();
          int max_s_waypoint = map_waypoints_s[map_waypoints_s.size() - 1];

          step += 1;

          int id_map_last = map_waypoints_x.size() - 1;
          int _close_way_point_id = ClosestWaypoint(car_x, car_y, map_waypoints_x, map_waypoints_y);

          int id_interp_start = _close_way_point_id - 5;
          int id_interp_end = _close_way_point_id + 8;

          vector<double> map_x_to_interp;
          vector<double> map_y_to_interp;
          vector<double> map_s_to_interp;

          double _map_s;
          double _map_x;
          double _map_y;

          for (int map_id = id_interp_start; map_id < id_interp_end; map_id++) {
            // 存疑
            if (map_id > id_map_last) {
              int _map_id = map_id - id_map_last - 1;
              _map_s = map_waypoints_s[_map_id] + max_s;
              _map_x = map_waypoints_x[_map_id];
              _map_y = map_waypoints_y[_map_id];
            } 
            // 存疑
            else if(map_id < 0) {
              int _map_id = id_map_last + map_id + 1;
              _map_s = map_waypoints_s[_map_id] - max_s;
              _map_x = map_waypoints_x[_map_id];              
              _map_y = map_waypoints_y[_map_id];              
            } else {
              _map_s = map_waypoints_s[map_id];
              _map_x = map_waypoints_x[map_id];
              _map_y = map_waypoints_y[map_id];
            }
            map_s_to_interp.push_back(_map_s);
            map_s_to_interp.push_back(_map_x);
            map_s_to_interp.push_back(_map_y);
          }

          tk::spline x_given_s;
          tk::spline y_given_s;
          x_given_s.set_points(map_s_to_interp, map_x_to_interp);
          y_given_s.set_points(map_s_to_interp, map_y_to_interp);

          vector<double> map_ss;
          vector<double> map_xs;
          vector<double> map_ys;

          double _s = map_s_to_interp[0];

          // 对地图曲线插值
          while(_s < map_s_to_interp[map_s_to_interp.size() - 1]) {
            double _x = x_given_s(_s);
            double _y = y_given_s(_s);
            map_ss.push_back(_s);
            map_xs.push_back(_x);
            map_ys.push_back(_y);
            _s += 0.1;
          }

          vector<Planner> planners;
          for (int i=0; i<3; i++) {
            double _target_d = 2.0 + 4* i - 0.15;
            Planner planner;
            MatrixXd s_trajectories(6, 0);
            VectorXd s_costs(0);
            MatrixXd d_trajectories(6, 0);
            VectorXd d_costs(0);

            planner.s_trajectories = s_trajectories;
            planner.s_costs = s_costs;
            planner.d_trajectories = d_trajectories;
            planner.d_costs = d_costs;
            planner.target_d = _target_d;
            planner.dist_to_target = 999.9;
            planner.obstacle_following = false;
            planner.feasible_traj_exist = true;
            planner.minimal_cost = 9999999.9;
            planner.optimal_s_id = 0;
            planner.optimal_d_id = 0;
            planner.iters = -1;
            planners.push_back(planner);
          }

          // 近距离范围内的他车
          vector<Vehicle> NearbyVehicles;
          for (int i=0; i<sensor_fusion.size(); i++) {
            // parsing the sensor fusion data
            // id, x, y, vx, vy, s, d
            double s_other = sensor_fusion[i][5];
            double s_dist = s_other - car_s;
            double d_other = sensor_fusion[i][6];
            // NEARBY VEHICLES
            double detect_range_front = 70.0;
            double detect_range_backward = 20.0;

            if ((s_dist < detect_range_front) && (s_dist >= - detect_range_backward) && (d_other > 0)) {
              Vehicle vehicle;
              vehicle.id = sensor_fusion[i][0];
              vehicle.x  = sensor_fusion[i][1];
              vehicle.y  = sensor_fusion[i][2];
              vehicle.vx = sensor_fusion[i][3];
              vehicle.vy = sensor_fusion[i][4];
              vehicle.s  = sensor_fusion[i][5];
              vehicle.d  = sensor_fusion[i][6];
              vehicle.speed = sqrt(vehicle.vx * vehicle.vx + vehicle.vy * vehicle.vy);

              NearbyVehicles.push_back(vehicle);
            }
          }

          int n_planning_horizon = 150;
          int n_pass = n_planning_horizon - prev_path_size;
          // int start_index = n_pass - 1;
          int start_index = n_pass - 1;
          double start_time = start_index * 0.02;

          double s0 = getPosition(optimal_s_coeff, start_time);
          if (s0 > max_s) {s0 = s0 - max_s;}
          // if (s0 > map_waypoints_s[map_waypoints_s.size()-1]) {s0 = s0 - max_s; cout << " [!!!!!] s0 = " << s0 << endl;}

          double s0dot = getVelocity(optimal_s_coeff, start_time);
          double s0ddot = getAcceleration(optimal_s_coeff, start_time);
          double d0 = getPosition(optimal_d_coeff, start_time);
          double d0dot = getVelocity(optimal_d_coeff, start_time);
          double d0ddot = getAcceleration(optimal_d_coeff, start_time);

          // cout << " [!!] s_offset = " << s_offset << endl;
          int mylane;
          // if (prev_path_size == 0) {mylane = getMyLane(car_d);}
          // else {mylane = getMyLane(d0);}
          mylane = getMyLane(car_d);

          for (int i=0; i<NearbyVehicles.size(); i++) {
            Vehicle _vehicle = NearbyVehicles[i];
            for (int j=0; j<planners.size(); j++) {
              // 当前车道内的他车
              if ((_vehicle.d > planners[j].target_d - 2) && (_vehicle.d <= planners[j].target_d + 2)) {

                // frontmost obstacle
                double from_ego_to_other = _vehicle.s - car_s;
                if (j == mylane) {
                  if (from_ego_to_other >= 3) {
                    planners[j].obstacles.push_back(_vehicle);
                  }
                }
                else {planners[j].obstacles.push_back(_vehicle);}
                // 保证target为距自车最近的他车
                if (from_ego_to_other >= -1.0){
                  if (from_ego_to_other < planners[j].dist_to_target) {
                    planners[j].dist_to_target = from_ego_to_other;
                    planners[j].target_to_follow = _vehicle;
                  }
                  if (from_ego_to_other <= 65) {
                    planners[j].obstacle_following = true;
                  }
                }
              }
            }
          }


          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }  // end "telemetry" if
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket if
  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  
  h.run();
}
