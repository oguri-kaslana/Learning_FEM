#include "fem_all.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

using namespace std;

double dot_vector(const vector<double>& a,
                  const vector<double>& b)
{
    const size_t n = min(a.size(), b.size());
    double s = 0.0;
    for (size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

double norm_vector(const vector<double>& a)
{
    return sqrt(max(0.0, dot_vector(a, a)));
}

void copy_vector(const vector<double>& src,
                 vector<double>& dst)
{
    dst = src;
}

void copy_matrix(const vector<double>& src,
                 vector<double>& dst)
{
    dst = src;
}

void zero_vector(vector<double>& x)
{
    fill(x.begin(), x.end(), 0.0);
}

void zero_matrix(vector<double>& A)
{
    fill(A.begin(), A.end(), 0.0);
}

void multiply_matrix_vector_dense(const vector<double>& A,
                                  const vector<double>& x,
                                  vector<double>& y,
                                  int n)
{
    y.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        const size_t base = static_cast<size_t>(i) * n;
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += A[base + j] * x[j];
        y[i] = s;
    }
}

void multiply_matrix_vector_sparse(const SparseCSR& A,
                                   const vector<double>& x,
                                   vector<double>& y)
{
    A.multiply_vector(x, y);
}

double calculate_residual_dense(const vector<double>& A,
                                const vector<double>& x,
                                const vector<double>& b,
                                int n)
{
    vector<double> Ax;
    multiply_matrix_vector_dense(A, x, Ax, n);
    double nr2 = 0.0;
    double nb2 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = Ax[i] - b[i];
        nr2 += r * r;
        nb2 += b[i] * b[i];
    }
    return sqrt(nr2) / max(1.0e-30, sqrt(nb2));
}

double calculate_residual_sparse(const SparseCSR& A,
                                 const vector<double>& x,
                                 const vector<double>& b)
{
    vector<double> Ax;
    A.multiply_vector(x, Ax);
    double nr2 = 0.0;
    double nb2 = 0.0;
    const int n = min(static_cast<int>(Ax.size()), static_cast<int>(b.size()));
    for (int i = 0; i < n; ++i) {
        const double r = Ax[i] - b[i];
        nr2 += r * r;
        nb2 += b[i] * b[i];
    }
    return sqrt(nr2) / max(1.0e-30, sqrt(nb2));
}

void swap_row(vector<double>& A,
              vector<double>& b,
              int n,
              int r1,
              int r2)
{
    if (r1 == r2) return;
    for (int j = 0; j < n; ++j) {
        swap(A[static_cast<size_t>(r1) * n + j], A[static_cast<size_t>(r2) * n + j]);
    }
    swap(b[r1], b[r2]);
}

double find_max_abs(const vector<double>& x)
{
    double m = 0.0;
    for (size_t i = 0; i < x.size(); ++i) m = max(m, fabs(x[i]));
    return m;
}

double calculate_max_displacement(const vector<double>& U)
{
    double m = 0.0;
    for (size_t i = 0; i + 1 < U.size(); i += 2) {
        const double umag = sqrt(U[i] * U[i] + U[i + 1] * U[i + 1]);
        if (umag > m) m = umag;
    }
    return m;
}

bool solve_gauss_dense(vector<double> A,
                       vector<double> b,
                       vector<double>& x,
                       int n)
{
    x.assign(n, 0.0);
    if (n <= 0 || static_cast<int>(A.size()) != n * n || static_cast<int>(b.size()) != n) return false;

    // 带部分主元高斯消元。
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double maxv = fabs(A[static_cast<size_t>(k) * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double v = fabs(A[static_cast<size_t>(i) * n + k]);
            if (v > maxv) {
                maxv = v;
                pivot = i;
            }
        }
        if (maxv < 1.0e-20) {
            cerr << "[solve_gauss_dense warning] singular matrix near row " << k << endl;
            return false;
        }
        swap_row(A, b, n, k, pivot);

        const double akk = A[static_cast<size_t>(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double factor = A[static_cast<size_t>(i) * n + k] / akk;
            A[static_cast<size_t>(i) * n + k] = 0.0;
            for (int j = k + 1; j < n; ++j) {
                A[static_cast<size_t>(i) * n + j] -= factor * A[static_cast<size_t>(k) * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) {
            s -= A[static_cast<size_t>(i) * n + j] * x[j];
        }
        const double aii = A[static_cast<size_t>(i) * n + i];
        if (fabs(aii) < 1.0e-20) return false;
        x[i] = s / aii;
    }

    return true;
}

bool solve_cholesky_dense_core(const vector<double>& A,
                               const vector<double>& b,
                               vector<double>& x,
                               int n)
{
    x.assign(n, 0.0);
    if (n <= 0 || static_cast<int>(A.size()) != n * n || static_cast<int>(b.size()) != n) return false;

    vector<double> L(static_cast<size_t>(n) * n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = A[static_cast<size_t>(i) * n + j];
            for (int k = 0; k < j; ++k) {
                s -= L[static_cast<size_t>(i) * n + k] * L[static_cast<size_t>(j) * n + k];
            }
            if (i == j) {
                if (s <= 1.0e-20) {
                    cerr << "[Cholesky warning] matrix is not positive definite at row " << i << endl;
                    return false;
                }
                L[static_cast<size_t>(i) * n + j] = sqrt(s);
            } else {
                L[static_cast<size_t>(i) * n + j] = s / L[static_cast<size_t>(j) * n + j];
            }
        }
    }

    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int k = 0; k < i; ++k) s -= L[static_cast<size_t>(i) * n + k] * y[k];
        y[i] = s / L[static_cast<size_t>(i) * n + i];
    }

    for (int i = n - 1; i >= 0; --i) {
        double s = y[i];
        for (int k = i + 1; k < n; ++k) s -= L[static_cast<size_t>(k) * n + i] * x[k];
        x[i] = s / L[static_cast<size_t>(i) * n + i];
    }

    return true;
}

void extract_free_dof_dense(const vector<double>& K,
                            const vector<double>& F,
                            const vector<int>& free_dof,
                            int n_dof,
                            vector<double>& Kff,
                            vector<double>& Ff)
{
    const int nf = static_cast<int>(free_dof.size());
    Kff.assign(static_cast<size_t>(nf) * nf, 0.0);
    Ff.assign(nf, 0.0);

    for (int i = 0; i < nf; ++i) {
        const int row = free_dof[i];
        Ff[i] = F[row];
        for (int j = 0; j < nf; ++j) {
            const int col = free_dof[j];
            Kff[static_cast<size_t>(i) * nf + j] = K[static_cast<size_t>(row) * n_dof + col];
        }
    }
}

void extract_free_dof_sparse(const SparseCSR& K,
                             const vector<double>& F,
                             const vector<int>& free_dof,
                             SparseCSR& Kff,
                             vector<double>& Ff)
{
    extract_free_csr(K, free_dof, Kff);
    extract_free_vector(F, free_dof, Ff);
}

void convert_sparse_to_dense_submatrix(const SparseCSR& K,
                                       const vector<int>& free_dof,
                                       vector<double>& Kff_dense)
{
    SparseCSR Kff;
    extract_free_csr(K, free_dof, Kff);
    convert_csr_dense(Kff, Kff_dense);
}

void extract_diag_dense(const vector<double>& A,
                        int n,
                        vector<double>& diag)
{
    diag.assign(n, 0.0);
    for (int i = 0; i < n; ++i) diag[i] = A[static_cast<size_t>(i) * n + i];
}

void extract_diag_sparse(const SparseCSR& A,
                         vector<double>& diag)
{
    diag.assign(A.n_row, 0.0);
    for (int i = 0; i < A.n_row; ++i) diag[i] = A.get_diag_value(i);
}
