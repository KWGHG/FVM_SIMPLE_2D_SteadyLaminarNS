#include "coefficient.h"

void buildMomentumCoeff(std::vector<T>& coefficients, Eigen::MatrixXd& Var_matrix, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::VectorXd& U_BCvector, Eigen::VectorXd& V_BCvector, Eigen::Vector2d& delta, double& nu, int& row, int& col, int& i, int& j)
{
	int indexI = 0, indexE = 0, indexN = 0, indexW = 0, indexS = 0, index_FUe = 0, index_FUw = 0, index_FVn = 0, index_FVs = 0;
	double AP = 0, APE = 0, APN = 0, APW = 0, APS = 0;
	double AE = 0, AW = 0, AN = 0, AS = 0;
	double deltaX, deltaY;
	Eigen::MatrixXi index(4, 2); index << i + 1, j, i, j + 1, i - 1, j, i, j - 1;
	deltaX = delta(0);
	deltaY = delta(1);

	indexI = (row - 1) * i + j;
	indexE = (row - 1) * index(0, 0) + index(0, 1);
	indexN = (row - 1) * index(1, 0) + index(1, 1);
	indexW = (row - 1) * index(2, 0) + index(2, 1);
	indexS = (row - 1) * index(3, 0) + index(3, 1);

	index_FUw = indexI;
	index_FUe = indexI + (row - 1);
	index_FVs = i + (col - 1) * j;
	index_FVn = index_FVs + (col - 1);

	if (i == 0) {
		if (j == 0) {//OK
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW =  nu * deltaY / (deltaX);//Velocity Inlet BC condition
			APS = nu * deltaX / (deltaY);//downWall BC condition
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;

			U_BCvector(indexI) = (0.01 * deltaY * 0.01 + nu * 0.01 * deltaY / (deltaX)) / Var_matrix.row(8)(indexI);//Velocity Inlet BC

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
		}
		else if (j == row - 2) {//ok
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = nu * deltaX / (deltaY);//upWall Bc condition
			APW = nu * deltaY / (deltaX);//Velocity Inlet BC condition
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			U_BCvector(indexI) = (0.01 * deltaY * 0.01 + nu * 0.01 * deltaY / (deltaX)) / Var_matrix.row(8)(indexI);//Velocity Inlet BC

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
		else {//ok
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW = nu * deltaY / (deltaX);//Velocity Inlet BC condition
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			U_BCvector(indexI) = (0.01 * deltaY * 0.01 + nu * 0.01 * deltaY / (deltaX)) / Var_matrix.row(8)(indexI);//Velocity Inlet BC

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
	}
	else if (i == col - 2) {
		//PressureOutlet no viscous effect
		if (j == 0) {
			APE = Var_matrix.row(0)(indexI) * deltaY;//pressureOutlet BC condition
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = nu * deltaX / (deltaY);//downWall BC condition
			AP = APE + APN + APW + APS;

			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;
			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
		}
		else if (j == row - 2) {
			APE = Var_matrix.row(0)(indexI) * deltaY;//pressureOutlet BC condition
			APN = nu * deltaX / (deltaY);//upWall BC condition
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
		else {
			APE = Var_matrix.row(0)(indexI) * deltaY;//pressureOutlet BC condition
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;
			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
	}
	else {
		if (j == 0) {
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = nu * deltaX / (deltaY);//downWall BC condition
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;
			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
		}
		else if (j == row - 2) {
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = nu * deltaX / (deltaY);//upWall BC condition
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
		else {
			APE = 0.5 * (fabs(Face_U.row(0)(index_FUe)) + Face_U.row(0)(index_FUe)) * deltaY + nu * deltaY / deltaX;
			APN = 0.5 * (fabs(Face_V.row(0)(index_FVn)) + Face_V.row(0)(index_FVn)) * deltaX + nu * deltaX / deltaY;
			APW = 0.5 * (fabs(Face_U.row(0)(index_FUw)) - Face_U.row(0)(index_FUw)) * deltaY + nu * deltaY / deltaX;
			APS = 0.5 * (fabs(Face_V.row(0)(index_FVs)) - Face_V.row(0)(index_FVs)) * deltaX + nu * deltaX / deltaY;
			AP = APE + APN + APW + APS;

			AE = -0.5 * (fabs(Face_U.row(0)(index_FUe)) - Face_U.row(0)(index_FUe)) * deltaY - nu * deltaY / deltaX;
			AN = -0.5 * (fabs(Face_V.row(0)(index_FVn)) - Face_V.row(0)(index_FVn)) * deltaX - nu * deltaX / deltaY;
			AW = -0.5 * (fabs(Face_U.row(0)(index_FUw)) + Face_U.row(0)(index_FUw)) * deltaY - nu * deltaY / deltaX;
			AS = -0.5 * (fabs(Face_V.row(0)(index_FVs)) + Face_V.row(0)(index_FVs)) * deltaX - nu * deltaX / deltaY;

			coefficients.emplace_back(T(indexI, indexI, AP / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexE, AE / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexN, AN / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexW, AW / Var_matrix.row(8)(indexI)));
			coefficients.emplace_back(T(indexI, indexS, AS / Var_matrix.row(8)(indexI)));
		}
	}
}

void PWIM(Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::Vector2d& delta, int& row, int& col, int& i, int& j)
{
	int indexI = 0, indexE = 0, indexEE = 0, indexN = 0, indexNN = 0, indexW = 0, indexS = 0, index_FUe = 0, index_FUw = 0, index_FVn = 0, index_FVs = 0;
	double deltaX, deltaY;
	deltaX = delta(0);
	deltaY = delta(1);

	Eigen::MatrixXi index(4, 2); index << i + 1, j, i, j + 1, i - 1, j, i, j - 1;

	indexI = (row - 1) * i + j;
	indexE = (row - 1) * index(0, 0) + index(0, 1);
	indexEE = (row - 1) * (index(0, 0) + 1) + index(0, 1);
	indexN = (row - 1) * index(1, 0) + index(1, 1);
	indexNN = (row - 1) * index(1, 0) + (index(1, 1) + 1);
	indexW = (row - 1) * index(2, 0) + index(2, 1);
	indexS = (row - 1) * index(3, 0) + index(3, 1);

	index_FUw = indexI;
	index_FUe = indexI + (row - 1);
	index_FVs = i + (col - 1) * j;
	index_FVn = index_FVs + (col - 1);

	if (i == 0) {
		if (j == 0) {
			Face_V.row(2)(index_FVs) = 0.0;
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 3) {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 2) {
			Face_V.row(2)(index_FVn) = 0.0;
		}
		else {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
			+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX / Var_matrix.row(8)(indexI)
				+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (Var_matrix.row(2)(indexEE) - Var_matrix.row(2)(indexI)) / deltaX / Var_matrix.row(8)(indexE)
				- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
		Face_U.row(2)(index_FUw) = 0.01;
	}
	else if (i == col - 3) {
		if (j == 0) {
			Face_V.row(2)(index_FVs) = 0.0;
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 3) {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 2) {
			Face_V.row(2)(index_FVn) = 0.0;
		}
		else {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
			+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexW)) / deltaX / Var_matrix.row(8)(indexI)
				+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (-Var_matrix.row(2)(indexW) - Var_matrix.row(2)(indexI)) / deltaX / Var_matrix.row(8)(indexE)
				- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
	}
	else if (i == col - 2) {
		if (j == 0) {
			Face_V.row(2)(index_FVs) = 0.0;
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 3) {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 2) {
			Face_V.row(2)(index_FVn) = 0.0;
		}
		else {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		Face_U.row(2)(index_FUe) = Var_matrix.row(6)(indexI);
	}
	else {
		if (j == 0) {
			Face_V.row(2)(index_FVs) = 0.0;
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 3) {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		else if (j == row - 2) {
			Face_V.row(2)(index_FVn) = 0.0;
		}
		else {
			Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
				+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY / Var_matrix.row(8)(indexI)
					+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY / Var_matrix.row(8)(indexN)
					- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
		}
		Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
			+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexW)) / deltaX / Var_matrix.row(8)(indexI)
				+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (Var_matrix.row(2)(indexEE) - Var_matrix.row(2)(indexI)) / deltaX / Var_matrix.row(8)(indexE)
				- (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
	}

	//if (i == 0) {//OK
	//	Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
	//		+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX
	//			+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (Var_matrix.row(2)(indexEE) - Var_matrix.row(2)(indexI)) / deltaX
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
	//}
	//else if (i == col - 3) { 
	//	Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
	//		+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexW)) / deltaX
	//			+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (-Var_matrix.row(2)(indexW) - Var_matrix.row(2)(indexI)) / deltaX
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
	//}
	//else if (i == col - 2) {
	//	Face_U.row(2)(index_FUe) = Var_matrix.row(6)(indexI);
	//}
	//else {
	//	Face_U.row(2)(index_FUe) = 0.5 * (Var_matrix.row(6)(indexI) + Var_matrix.row(6)(indexE))
	//		+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexW)) / deltaX
	//			+ 0.5 * LHSD_inverse.coeffRef(indexE, indexE) * (Var_matrix.row(2)(indexEE) - Var_matrix.row(2)(indexI)) / deltaX
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE)) * (Var_matrix.row(2)(indexE) - Var_matrix.row(2)(indexI)) / deltaX);
	//}

	//if (j == 0) {
	//	Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
	//		+ 0.5 * deltaX * deltaY *(0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY
	//			+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexN, indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
	//}
	//else if(j == row - 3){
	//	Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
	//		+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY
	//			+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexS)) / deltaY
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexN, indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
	//}
	//else if (j == row - 2) {
	//	Face_V.row(2)(index_FVn) = 0.0;
	//}
	//else {
	//	Face_V.row(2)(index_FVn) = 0.5 * (Var_matrix.row(7)(indexI) + Var_matrix.row(7)(indexN))
	//		+ 0.5 * deltaX * deltaY * (0.5 * LHSD_inverse.coeffRef(indexI, indexI) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexS)) / deltaY
	//			+ 0.5 * LHSD_inverse.coeffRef(indexN, indexN) * (Var_matrix.row(2)(indexNN) - Var_matrix.row(2)(indexI)) / deltaY
	//			- (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexN, indexN)) * (Var_matrix.row(2)(indexN) - Var_matrix.row(2)(indexI)) / deltaY);
	//}
}

void buildPoissonCoeff(std::vector<T>& coefficients, Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::VectorXd& RHSP, Eigen::Vector2d& delta, int& row, int& col, int& i, int& j)
{
	int indexI = 0, indexE = 0, indexN = 0, indexW = 0, indexS = 0, index_FUe = 0, index_FUw = 0, index_FVn = 0, index_FVs = 0;
	double AP = 0, APE = 0, APN = 0, APW = 0, APS = 0;
	double AE = 0, AW = 0, AN = 0, AS = 0;
	double source = 0;
	double deltaX, deltaY;
	Eigen::MatrixXi index(4, 2); index << i + 1, j, i, j + 1, i - 1, j, i, j - 1;
	deltaX = delta(0);
	deltaY = delta(1);

	indexI = (row - 1) * i + j;
	indexE = (row - 1) * index(0, 0) + index(0, 1);
	indexN = (row - 1) * index(1, 0) + index(1, 1);
	indexW = (row - 1) * index(2, 0) + index(2, 1);
	indexS = (row - 1) * index(3, 0) + index(3, 1);

	index_FUw = indexI;
	index_FUe = indexI + (row - 1);
	index_FVs = i + (col - 1) * j;
	index_FVn = index_FVs + (col - 1);

	//APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE));
	//APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexN, indexN));
	//APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexW, indexW));
	//APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexS, indexS));

	//AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE));
	//AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexN, indexN));
	//AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexW, indexW));
	//AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexS, indexS));

	source = -((Face_U.row(2)(index_FUe) - Face_U.row(2)(index_FUw)) * deltaY + (Face_V.row(2)(index_FVn) - Face_V.row(2)(index_FVs)) * deltaX);
	RHSP(indexI) = source;

	if (i == 0) {
		if (j == 0) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));

			AP = APE + APN;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexN, AN));
		}
		else if (j == row - 2) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
		else {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APN + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexN, AN));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
	}
	else if (i == col - 2) {
		if (j == 0) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));//pressureOutlet
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));

			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW))
				+ 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));

			AP = APE + APN + APW;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexN, AN));
			coefficients.emplace_back(T(indexI, indexW, AW));
		}
		else if (j == row - 2) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));//pressureOutlet
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW))
				+ 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APW + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexW, AW));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
		else {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));//pressureOutlet
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW))
				+ 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APN + APW + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexN, AN));
			coefficients.emplace_back(T(indexI, indexW, AW));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
	}
	else {
		if (j == 0) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));

			AP = APE + APN + APW;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexN, AN));
			coefficients.emplace_back(T(indexI, indexW, AW));
		}
		else if (j == row - 2) {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APW + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexW, AW));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
		else {
			APE = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			APN = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			APW = 0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			APS = 0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AE = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE));
			AN = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN));
			AW = -0.5 * pow(deltaY, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexW, indexW) / Var_matrix.row(8)(indexW));
			AS = -0.5 * pow(deltaX, 2) * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexS, indexS) / Var_matrix.row(8)(indexS));

			AP = APE + APN + APW + APS;
			coefficients.emplace_back(T(indexI, indexI, AP));
			coefficients.emplace_back(T(indexI, indexE, AE));
			coefficients.emplace_back(T(indexI, indexN, AN));
			coefficients.emplace_back(T(indexI, indexW, AW));
			coefficients.emplace_back(T(indexI, indexS, AS));
		}
	}
}

void gradientP(Eigen::VectorXd& RHSU, Eigen::VectorXd& RHSV, Eigen::MatrixXd& Var_matrix, Eigen::SparseMatrix<double>& LHSD_inverse, Eigen::MatrixXd& Face_U, Eigen::MatrixXd& Face_V, Eigen::Vector2d& delta, double& nu, double& alpha_velocity, int& row, int& col, int& i, int& j, int state)
{
	int indexI = 0, indexE = 0, indexN = 0, indexW = 0, indexS = 0, index_FUe = 0, index_FUw = 0, index_FVn = 0, index_FVs = 0;
	double deltaX, deltaY;
	Eigen::MatrixXi index(4, 2); index << i + 1, j, i, j + 1, i - 1, j, i, j - 1;
	deltaX = delta(0);
	deltaY = delta(1);

	indexI = (row - 1) * i + j;
	indexE = (row - 1) * index(0, 0) + index(0, 1);
	indexN = (row - 1) * index(1, 0) + index(1, 1);
	indexW = (row - 1) * index(2, 0) + index(2, 1);
	indexS = (row - 1) * index(3, 0) + index(3, 1);

	index_FUw = indexI;
	index_FUe = indexI + (row - 1);
	index_FVs = i + (col - 1) * j;
	index_FVn = index_FVs + (col - 1);

	if (state == 0) {
		if (i == 0) {
			if (j == 0) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else if (j == row - 2) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexI)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
		}
		else if (i == col - 2) {
			if (j == 0) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) + Var_matrix.row(2)(indexI)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else if (j == row - 2) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) + Var_matrix.row(2)(indexI)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexI)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) + Var_matrix.row(2)(indexI)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
		}
		else {
			if (j == 0) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexI) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else if (j == row - 2) {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexI)) * deltaX / Var_matrix.row(8)(indexI);
			}
			else {
				RHSU(indexI) = 0.5 * (Var_matrix.row(2)(indexW) - Var_matrix.row(2)(indexE)) * deltaY / Var_matrix.row(8)(indexI);
				RHSV(indexI) = 0.5 * (Var_matrix.row(2)(indexS) - Var_matrix.row(2)(indexN)) * deltaX / Var_matrix.row(8)(indexI);
			}
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////
	if (state == 1) {
		if (i == 0) {
			if (j == 0) {
				Face_V.row(1)(index_FVs) = 0.0;
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn) 
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;
						
			}
			else if (j == row - 2) {
				Face_V.row(1)(index_FVn) = 0.0;
			}
			else {
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn)
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;
			}
			Face_U.row(1)(index_FUe) = Face_U.row(2)(index_FUe)
				+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI / Var_matrix.row(8)(indexI)) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE)) * (Var_matrix.row(5)(indexE) - Var_matrix.row(5)(indexI)) * deltaY;
			Face_U.row(1)(index_FUw) = 0.01;
		}
		else if (i == col - 2) {
			if (j == 0) {
				Face_V.row(1)(index_FVs) = 0.0;
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn)
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;

			}
			else if (j == row - 2) {
				Face_V.row(1)(index_FVn) = 0.0;
			}
			else {
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn)
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;
			}
			Face_U.row(1)(index_FUe) = Var_matrix.row(6)(indexI);
		}
		else {
			if (j == 0) {
				Face_V.row(1)(index_FVs) = 0.0;
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn)
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;

			}
			else if (j == row - 2) {
				Face_V.row(1)(index_FVn) = 0.0;
			}
			else {
				Face_V.row(1)(index_FVn) = Face_V.row(2)(index_FVn)
					+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexN, indexN) / Var_matrix.row(8)(indexN)) * (Var_matrix.row(5)(indexN) - Var_matrix.row(5)(indexI)) * deltaX;
			}
			Face_U.row(1)(index_FUe) = Face_U.row(2)(index_FUe)
				+ alpha_velocity * 0.5 * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI) + LHSD_inverse.coeffRef(indexE, indexE) / Var_matrix.row(8)(indexE)) * (Var_matrix.row(5)(indexE) - Var_matrix.row(5)(indexI)) * deltaY;
			//std::cout << "FU " << Face_U.row(2)(index_FUe) << "\t" << (LHSD_inverse.coeffRef(indexI, indexI) + LHSD_inverse.coeffRef(indexE, indexE)) << "\t" << (Var_matrix.row(5)(indexE) - Var_matrix.row(5)(indexI)) << "\n";
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	if(state == 2){
		if (i == 0) {
			if (j == 0) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else if (j == row - 2) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexI)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
		}
		else if (i == col - 2) {
			if (j == 0) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) + Var_matrix.row(5)(indexI)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else if (j == row - 2) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) + Var_matrix.row(5)(indexI)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexI)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) + Var_matrix.row(5)(indexI)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
		}
		else {
			if (j == 0) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexI) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else if (j == row - 2) {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexI)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			else {
				Var_matrix.row(3)(indexI) = Var_matrix.row(6)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexW) - Var_matrix.row(5)(indexE)) * deltaY * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
				Var_matrix.row(4)(indexI) = Var_matrix.row(7)(indexI) + alpha_velocity * 0.5 * (Var_matrix.row(5)(indexS) - Var_matrix.row(5)(indexN)) * deltaX * (LHSD_inverse.coeffRef(indexI, indexI) / Var_matrix.row(8)(indexI));
			}
			
		}
	}
}