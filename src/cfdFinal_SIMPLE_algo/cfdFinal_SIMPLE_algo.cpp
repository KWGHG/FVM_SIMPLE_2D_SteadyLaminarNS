#include "coefficient.h"
int main()
{
    int maxIteration = 100;

    //i dir is left to rignt, j dir is lower to upper
    double rho = 1.205, nu = 1.82E-5, width = 0.01, length = 0.1, alpha_pressure = 0.01, alpha_velocity = 0.5;
    int row = 0, col = 0;
    std::string filename = "";
    std::cout << "Input filename (only filename!!)" << std::endl;
    std::cin >> filename;
    std::cout << "Input maxIteration, alpha_pressure, alpha_velocity" << std::endl;
    std::cin >> maxIteration >> alpha_pressure >> alpha_velocity;
    std::ifstream ifs(filename + ".dat");
    ifs >> row >> col;
    int gridSizes = (row - 1) * (col - 1);
    std::cout << row << "\t" << col << "\n";
    Eigen::Tensor<double, 3> grid(2, col - 1, row - 1); // Eigen::Tensor<int, 3> A(a, b, c)a=變數, 0:x網格值, 1: y網格值, 2: Volume; b=x方向網格數; c=y方向網格數
    Eigen::Tensor<double, 3> dots(2, col, row); // Eigen::Tensor<int, 3> A(a, b, c)a=變數, 0: x網格值, 1: y網格值; b=x方向網格數; c=y方向網格數
    grid.setConstant(0);
    dots.setConstant(0);

    for (int i = 0; i < col; ++i) {
        for (int j = 0; j < row; ++j) {
            ifs >> dots(0, i, j) >> dots(1, i, j);
            //std::cout << "(" << i << ", " << j << ") " << dots(0, i, j) << "\t" << dots(1, i, j) << std::endl;
        }
    }
    ifs.close();

    Eigen::MatrixXd Var_matrix = Eigen::MatrixXd::Zero(9, gridSizes);// Unow, Vnow, Pnow, Unew, Vnew, Pnew, Upredict, Vpredict, Volume
    int indexI = 0, indexJ = 0;
    for (int i = 0; i < col - 1; ++i) {
        for (int j = 0; j < row - 1; ++j) {
            indexI = (row - 1) * i + j;
            grid(0, i, j) = 0.25 * (dots(0, i, j) + dots(0, i + 1, j) + dots(0, i, j + 1) + dots(0, i + 1, j + 1));
            grid(1, i, j) = 0.25 * (dots(1, i, j) + dots(1, i + 1, j) + dots(1, i, j + 1) + dots(1, i + 1, j + 1));
            Var_matrix.row(8)(indexI) = 1 / ((dots(0, i + 1, j) - dots(0, i, j)) * (dots(1, i, j + 1) - dots(1, i, j)));
        }
    }

    //B.C conditions
    Eigen::Vector3d velocityInlet(0.01, grid(0, 0, 0), 2.0);// BC Velocity inlet (velocity, distance between face to cell center, west face)
    Eigen::Vector3d downWall(0, grid(1, 0, 0), 3.0);//BC down Wall (velocity, distance, south face)
    Eigen::Vector3d upWall(0, width - grid(1, col - 2, row - 2), 1.0);//BC up Wall (velocity, distance, north face)
    Eigen::Vector3d pressureOutlet(0, length - grid(0, col - 2, row - 2), 0.0);//BC pressure outlet (Pressure, distance between face to cell center, east face)

    //Variables
    
    Eigen::MatrixXd Face_U = Eigen::MatrixXd::Zero(3, gridSizes + (row - 1));// Un, Unew, Upredict
    Eigen::MatrixXd Face_V = Eigen::MatrixXd::Zero(3, gridSizes + (col - 1));// Vn, Vnew, Vpredict
    Eigen::MatrixXd Residual = Eigen::MatrixXd::Zero(3, gridSizes);//Residual matrix
    Eigen::SparseMatrix<double> LHS(gridSizes, gridSizes);//Momentum matrix(LHS)
    Eigen::SparseMatrix<double> LHSD(gridSizes, gridSizes);//LHS diagonal matrix
    Eigen::SparseMatrix<double> LHSD_inverse(gridSizes, gridSizes);//inverse LHS diagonal matrix
    Eigen::SparseMatrix<double> LHSP(gridSizes, gridSizes);//Poisson matrix(LHS)


    Eigen::VectorXd RHSU = Eigen::VectorXd::Zero(gridSizes);// X-Momentum RHS
    Eigen::VectorXd RHSV = Eigen::VectorXd::Zero(gridSizes);// Y-Momentum RHS
    Eigen::VectorXd RHSP = Eigen::VectorXd::Zero(gridSizes);// Poisson RHS
    Eigen::VectorXd U_bc = Eigen::VectorXd::Zero(gridSizes);// return U_BC  condition to add RHS
    Eigen::VectorXd V_bc = Eigen::VectorXd::Zero(gridSizes);// return V_BC  condition to add RHS

    std::vector<T> LHS_coefficients;
    std::vector<T> LHSP_coefficients;

    Eigen::Vector2d delta(grid(0, 1, 0) - grid(0, 0, 0), grid(1, 0, 1) - grid(1, 0, 0));

    //Matrix solver
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
    //cg.setMaxIterations(10000);
    //Initial conditions
    for (int i = 0; i < gridSizes; ++i) {
        Var_matrix.row(0)(i) = velocityInlet(0);
    }
    for (int i = 0; i < Face_U.cols(); ++i) {
        if (i < row - 1)
        {
            Face_U.row(0)(i) = velocityInlet(0);
            Face_U.row(1)(i) = velocityInlet(0);
            Face_U.row(2)(i) = velocityInlet(0);
        }
        else {
            Face_U.row(0)(i) = velocityInlet(0);
        }
    }

    double residual_maxU = 0, residual_maxV = 0, residual_maxP = 0;
    double first_residual_maxU = 1, first_residual_maxV = 1, first_residual_maxP = 1;

    Eigen::SparseMatrix<double> LHS_(gridSizes, gridSizes);//Momentum matrix(LHS)
    Eigen::VectorXd RHSU_ = Eigen::VectorXd::Zero(gridSizes);// X-Momentum RHS
    Eigen::VectorXd RHSV_ = Eigen::VectorXd::Zero(gridSizes);// Y-Momentum RHS
    ////////////////////////////////////////////////////////////////////////////////////////////
    for(int iter = 0 ; iter < maxIteration ; ++iter){
        //Solving momentum equation
        LHS_coefficients.clear();
        LHSD_inverse.setZero();
        LHS.setZero();

        for (int i = 0; i < col - 1; ++i) {
            for (int j = 0; j < row - 1; ++j) {
                buildMomentumCoeff(LHS_coefficients, Var_matrix, Face_U, Face_V, U_bc, V_bc, delta, nu, row, col, i, j);
                gradientP(RHSU, RHSV, Var_matrix, LHSD_inverse, Face_U, Face_V, delta, nu, alpha_velocity, row, col, i, j, 0);
            }
        }

        LHS.setFromTriplets(LHS_coefficients.begin(), LHS_coefficients.end());
        RHSU = RHSU + U_bc;//Corrector BC conditions
        RHSV = RHSV + V_bc;//Corrector BC conditions

        
        cg.compute(LHS);//Solving predicted velocity
        Var_matrix.row(6).transpose() = cg.solve(RHSU);
        Var_matrix.row(7).transpose() = cg.solve(RHSV);

        //std::cout << "\n==========P_Velocity==========\n" << "\n";
        //for (int i = 0; i < gridSizes; ++i) {
        //    std::cout << Var_matrix.row(6)(i) << "\t" << Var_matrix.row(7)(i) << "\n";
        //}

        //Get inverse LHS diagonal matrix
        LHSD.setZero();
        LHSD_inverse.setZero();
        LHSD = LHS.diagonal().asDiagonal();
        for (int i = 0; i < gridSizes; ++i) {
            LHSD_inverse.coeffRef(i, i) = 1 / LHSD.coeffRef(i, i);
        }

        //Calculating predicted face velocity that use PWIM
        for (int i = 0; i < col - 1; ++i) {
            for (int j = 0; j < row - 1; ++j) {
                PWIM(Var_matrix, LHSD_inverse, Face_U, Face_V, delta, row, col, i, j);
            }
        }

        //Calculate Poisson equation to get modify pressure
        LHSP_coefficients.clear();
        for (int i = 0; i < col - 1; ++i) {
            for (int j = 0; j < row - 1; ++j) {
                buildPoissonCoeff(LHSP_coefficients, Var_matrix, LHSD_inverse, Face_U, Face_V, RHSP, delta, row, col, i, j);
            }
        }
        LHSP.setFromTriplets(LHSP_coefficients.begin(), LHSP_coefficients.end());

        cg.compute(LHSP);
        Var_matrix.row(5).transpose() = cg.solve(RHSP);

        Var_matrix.row(5) = alpha_pressure * Var_matrix.row(5);
        for (int i = 0; i < col - 1; ++i) {
            for (int j = 0; j < row - 1; ++j) {
                gradientP(RHSU, RHSV, Var_matrix, LHSD_inverse, Face_U, Face_V, delta, nu, alpha_velocity, row, col, i, j, 2);
                gradientP(RHSU, RHSV, Var_matrix, LHSD_inverse, Face_U, Face_V, delta, nu, alpha_velocity, row, col, i, j, 1);
            }
        }

        Residual.row(0).transpose() = RHSU - LHS * Var_matrix.row(3).transpose();
        Residual.row(1).transpose() = RHSV - LHS * Var_matrix.row(4).transpose();
        Residual.row(2).transpose() = RHSP - LHSP * Var_matrix.row(5).transpose();

        residual_maxU = Residual.row(0).maxCoeff();
        residual_maxV = Residual.row(1).maxCoeff();
        residual_maxP = Residual.row(2).maxCoeff();

        if (iter == 0) {
            first_residual_maxU = residual_maxU;
            first_residual_maxV = residual_maxV;
            first_residual_maxP = residual_maxP;
        }
        else {
            residual_maxU = residual_maxU / first_residual_maxU;
            residual_maxV = residual_maxV / first_residual_maxV;
            residual_maxP = residual_maxP / first_residual_maxP;
        }

        std::cout << "\n" << iter <<"==================Residual================\n";
        std::cout << "residual_maxU: " << residual_maxU << "\n";
        std::cout << "\nresidual_maxV: " << residual_maxV << "\n";
        std::cout << "\nresidual_maxP: " << residual_maxP << "\n";

        Var_matrix.row(5) = Var_matrix.row(2) + Var_matrix.row(5);
        Var_matrix.row(0) = Var_matrix.row(3);
        Var_matrix.row(1) = Var_matrix.row(4);
        Var_matrix.row(2) = Var_matrix.row(5);
        Face_U.row(0) = Face_U.row(1);
        Face_V.row(0) = Face_V.row(1);
      
}
    
    std::ofstream ofs1(filename + ".oup");
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(col, row);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(col, row);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(col, row);
    ofs1 << "X [m]" << "\t" << "Y [m]" << "\t" << "u[m/s]" << "\t" << "v [m/s]" << "\t" << "P [Pa]\n";
    int A, B, C, D;
    for (int i = 0; i < col-1; ++i) {
        for (int j = 0; j < row-1 ; ++j) {
            if (i < col - 2) {
                if (j < row - 2) {
                    A = (row - 1) * i + j;
                    B = (row - 1) * (i + 1) + j;
                    C = (row - 1) * i + (j + 1);
                    D = (row - 1) * (i + 1) + (j + 1);
                    U(i + 1, j + 1) = 0.25 * (Var_matrix.row(0)(A) + Var_matrix.row(0)(B) + Var_matrix.row(0)(C) + Var_matrix.row(0)(D));
                    V(i + 1, j + 1) = 0.25 * (Var_matrix.row(1)(A) + Var_matrix.row(1)(B) + Var_matrix.row(1)(C) + Var_matrix.row(1)(D));
                    P(i + 1, j + 1) = 0.25 * (Var_matrix.row(2)(A) + Var_matrix.row(2)(B) + Var_matrix.row(2)(C) + Var_matrix.row(2)(D));
                }  
            }
        }
    }

    for (int i = 0; i < col; ++i) {
        for (int j = 0; j < row; ++j) {
            P(0, 0) = 1;
            if (i == 0) {
                if (j == 0) {
                    P(0, 0) =P(1, 1);
                }
                else if (j == row - 1) {
                    P(0, j) = P(1, j - 1);
                }
                else {
                    U(i, j) = 0.01;
                    V(i, j) = 0.0;
                    P(i, j) = P(i + 1, j);
                }
            }
            else if (i == col - 1 && j != 0 && j != row - 1) {
                if (j == 0) {
                    P(i, j) = P(i, j + 1);
                }
                else if (j == row - 1) {
                    P(i, j) = P(i, j - 1);
                }
                else {
                    U(i, j) = U(i - 1, j);
                    V(i, j) = V(i - 1, j);
                }
            }
            else {
                if (j == 0) {
                    P(i, j) = P(i, j + 1);
                }
                if (j == row - 1) {
                    P(i, j) = P(i, j - 1);
                }
            }
            
            ofs1 << dots(0, i, j) << "\t" << dots(1, i, j) << "\t" << U(i, j) << "\t" << V(i, j) << "\t" << P(i, j) << "\n";
        }
    }
    
    
    ofs1.close();
    system("pause");
    return 0;
}
