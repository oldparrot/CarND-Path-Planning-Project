#pragma once

#include <vector>
#include <map>
#include <string>

class Vehicle {
 public:
    double s;
    double s_d;
    double s_dd;
    double d;
    double d_d;
    double d_dd;
    std::string state;
    std::vector<std::string> available_states;
    std::vector<double> s_traj_coeffs;
    std::vector<double> d_traj_coeffs;

    Vehicle();
    Vehicle(double s, double s_d, double s_dd, double d, double d_d, double d_dd);

    virtual ~Vehicle();
};