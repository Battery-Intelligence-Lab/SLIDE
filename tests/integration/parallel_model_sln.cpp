#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>

using namespace std; // bu yerine fonksiyonun içine ekle
using namespace Eigen;

double ocv_eval(const vector<double>& ocv_coefs, double z) {
    double ocv = 0.0;
    for (size_t i = 0; i < ocv_coefs.size(); ++i) {
        ocv += ocv_coefs[i] * pow(z, ocv_coefs.size() - 1 - i);
    }
    return ocv;
} // already existing function --> otomatik olarak değer fonksiyona gelecek

void new_compute_A11_A12_A21_A22(int n_par, VectorXd R, VectorXd C, VectorXd Q, VectorXd tau, VectorXd r, MatrixXd &A11, MatrixXd &A12, MatrixXd &A21, MatrixXd &A22, MatrixXd &m) {
    // Initialize A11 and A12
    A11 = MatrixXd::Zero(2 * n_par, 2 * n_par);
    A12 = MatrixXd::Zero(2 * n_par, 2 * n_par);
    
    A11.block(0, 0, 2, 2) << 0, 0, 0, -tau(0);
    A12.block(0, 0, 2, 2) << 1/Q(0), 0, 0, 1/C(0);

    for (int i = 1; i < n_par; ++i) {
        A11.block(2 * i, 2 * i, 2, 2) << 0, 0, 0, -tau(i);
        A12.block(2 * i, 2 * i, 2, 2) << 1/Q(i), 0, 0, 1/C(i);
    }

    // Initialize A22
    A22 = MatrixXd::Zero(n_par, n_par);
    VectorXd main_diag = -(r + R);
    main_diag(0) = r(0);

    VectorXd lower_diag = r.head(n_par - 1);
    lower_diag(0) = r(0);

    A22.diagonal() = main_diag;
    A22.diagonal(-1) = lower_diag;

    for (int i = 1; i < n_par; ++i) {
        for (int j = i + 1; j < n_par; ++j) {
            A22(i, j) = -R(i);
        }
    }

    A22.row(0).setOnes();
    cout << "A22:\n" << A22 << endl;

    // Initialize A21
    A21 = MatrixXd::Zero((n_par - 1) * 2, n_par);
    for (int i = 0; i < n_par - 1; ++i) {
        A21.block(2 * i, 0, 2, n_par) = MatrixXd::Identity(2, n_par);
    }

    // Calculate m (inverse of A22)
    m = MatrixXd::Zero(n_par, n_par);
    double Rsum_inv = 1 / R.sum();

    for (int ktil = 0; ktil < n_par; ++ktil) {
        for (int i = 0; i < n_par - 1; ++i) {
            if (ktil - i - 1 == 0) {
                m(i + 1, i) = pow(1 / R(i + 1), 2) * Rsum_inv - 1 / R(i + 1);
            } else {
                m(ktil, i) = 1 / (R(ktil) * R(i + 1)) * Rsum_inv;
            }
        }
    }

    m(n_par - 1, n_par - 1) = 1 / (R(n_par - 1) * Rsum_inv);
    for (int i = 0; i < n_par - 1; ++i) {
        m(i, n_par - 1) = 1 / (R(i) * Rsum_inv);
    }

    MatrixXd inv_A22 = A22.inverse();
    double error = (inv_A22 - m).norm();

    cout << "Error: " << error << endl;
}

int main() {
    int n_par = 3;
    VectorXd R(n_par), C(n_par), Q(n_par), tau(n_par), r(n_par);

    R << 0, 0, 0;
    C << 634.0, 634.0, 634.0;
    Q << 9000.0, 9000.0, 9000.0;
    tau << 0.04, 0.04, 0.04;
    r << 0.029, 0.029, 0.029;

    MatrixXd A11, A12, A21, A22, m;
    new_compute_A11_A12_A21_A22(n_par, R, C, Q, tau, r, A11, A12, A21, A22, m);

    cout << "A11:\n" << A11 << endl;
    cout << "A12:\n" << A12 << endl;
    cout << "A21:\n" << A21 << endl;
    cout << "A22:\n" << A22 << endl;
    cout << "m:\n" << m << endl;

    return 0;
}
