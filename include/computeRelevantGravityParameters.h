/// @Author       : linfd 3039562364@qq.com
/// @Date         : 2024-10-21 08:02:37
/// @LastEditTime : 2024-10-22 11:19:18
/// @FilePath     : \computeRelevantGravityParameters\include\computeRelevantGravityParameters.h
/// @Description  :

#pragma once
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <Python.h>
#include <vector>
#include <format>
using namespace std;

const double PI = 3.14159265358979;
#define DEG2RAD(x) ((x) * (PI) / 180.0)

tuple<double, double, double> BLH2XYZ(double B, double L, double H, double a, double f);
tuple<double, double, double> XYZ2R_THETA_LAMBDA(double X, double Y, double Z);

class ReferenceEllipsoid
{
public:
    ReferenceEllipsoid(const string name, double a, double f, double GM, double omega);
    double getA() { return m_a; }
    double getF() { return m_f; }
    void outputReferenceEllipsoidParameters();
    bool readCmn_AndSmn_(const string &fileName);
    bool computePmn_(double theta, vector<vector<double>> &Pnm_);
    double computeNormalGravityOnEllipsoidalSurface(double phi);
    double computeNormalGravityUpEllipsoidalSurface(double phi, double h);
    double computeGravityDisturbance(double theta, double lambda, double r, vector<vector<double>> Pnm_);
    double computeGravityAnomaly(double theta, double lambda, double r, vector<vector<double>> Pnm_);
    double computePerturbationPosition(double theta, double lambda, double r, vector<vector<double>> Pnm_);
    double computeGeoidDifference(double theta, double lambda, double r, vector<vector<double>> Pnm_);
    tuple<double, double, double, double> computeGeodeticGravity(double theta, double lambda, double r, vector<vector<double>> Pnm_);

private:
    string m_name;
    double m_a;
    double m_f;
    double m_b;
    double m_GM;
    double m_omega;
    double m_J2_;
    double m_e2;
    double m_ee2;
    double m_E;
    double m_c;
    double m_U0;
    double m_m;
    double m_garmaA;
    double m_garmaB;
    int m_order = 360;
    vector<vector<double>> m_Cnm_;
    vector<vector<double>> m_Snm_;
};