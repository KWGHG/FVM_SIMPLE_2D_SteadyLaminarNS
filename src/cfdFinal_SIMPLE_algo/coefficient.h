#pragma once
#include <iostream>
#include <fstream>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Sparse>
#include<Eigen/IterativeLinearSolvers>	
#include <vector>
#include <cmath>
typedef Eigen::Triplet<double> T;

//把邊界的值存入BC_Vector裡，BC_Vector(大小為所有邊界網格的數量)再帶回到gradientP。
void buildMomentumCoeff(std::vector<T>& coefficients, Eigen::MatrixXd& Var_matrix, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::VectorXd& U_BCvector, Eigen::VectorXd& V_BCvector, Eigen::Vector2d& delta, double& nu, int& row, int& col, int& i, int& j);
void PWIM(Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::Vector2d& delta, int& row, int& col, int& i, int& j);
void buildPoissonCoeff(std::vector<T>& coefficients, Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::VectorXd& RHSP, Eigen::Vector2d& delta, int& row, int& col, int& i, int& j);
void gradientP(Eigen::VectorXd& RHSU, Eigen::VectorXd& RHSV, Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::Vector2d& delta, double& nu, double& alpha_velocity, int& row, int& col, int& i, int& j, int state);