#include "lrs.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "score.h"

namespace colors {

inline int count_larger_than(const Vector& v, double value) {
  int count = 0;
  for (int i = 0; i < v.size(); ++i) {
    if (v[i] > value) ++count;
  }
  return count;
}

int calc_lrs(const double lambda, const Matrix& M, Matrix& L, Matrix& S) {
  int dim = M.rows();

  Matrix Y = M;  // copy
  L = Matrix::Zero(dim, dim);
  S = Matrix::Zero(dim, dim);
  Array Z = Array::Zero(dim, dim);

  Eigen::JacobiSVD<Matrix> svd_only_singlar_values(Y);
  const double norm_two = svd_only_singlar_values.singularValues()(0);
  const double norm_inf = Y.array().abs().maxCoeff() / lambda;
  const double dual_norm = std::max(norm_two, norm_inf);
  const double d_norm = M.norm();

  Y /= dual_norm;

  double mu = 1.25 / norm_two;
  const double rho = 1.5;
  const double mu_bar = mu * 1.0e+7;

  bool converged = false;
  int max_iter = 1000;
  double error_tolerance = 1.0e-7;
  int total_svd = 0;
  int sv = 10;
  int iter = 0;

  while (++iter < max_iter) {
    // update sparse matrix S
    Array temp_T = M - L + (1.0 / mu) * Y;
    S = (temp_T - lambda / mu).max(Z) + (temp_T + lambda / mu).min(Z);
    S = S.array().max(Z);

    // update low-rank matrix L
    Eigen::JacobiSVD<Matrix> svd(M - S + 1.0 / mu * Y,
                                 Eigen::ComputeFullU | Eigen::ComputeFullV);
    Matrix U = svd.matrixU();
    Matrix V = svd.matrixV();
    Vector singular_values = svd.singularValues();

    int svp = count_larger_than(singular_values, 1.0 / mu);
    if (svp < sv) {
      sv = std::min(svp + 1, dim);
    } else {
      sv = std::min(svp + static_cast<int>(0.05 * dim + 0.5), dim);
    }

    Matrix S_th =
        (singular_values.head(svp).array() - 1.0 / mu).matrix().asDiagonal();
    L = U.leftCols(svp) * S_th * V.leftCols(svp).transpose();
    L = L.array().max(Z);

    total_svd += 1;
    Matrix D = M - L - S;
    Y = Y + mu * D;

    // update mu
    mu = std::min(mu * rho, mu_bar);

    double objective = D.norm() / d_norm;
    if (objective < error_tolerance) {
      break;
    }
  }
  return 0;
}

int calc_lrs(const double lambda, const Matrix& M, const char* output_path) {
  Matrix L, S;
  calc_lrs(lambda, M, L, S);

  Matrix apc(M.cols(), M.cols());
  calc_apc(M, apc);

  save_lrs(M, S, apc, output_path);

  return 0;
}

int save_lrs(const Matrix& M, const Matrix& S, const Matrix& apc,
             const char* output_path) {
  std::ofstream fout;
  fout.open(output_path);
  if (fout.fail()) {
    std::cerr << "Failed in open file " << output_path << std::endl;
  }

  fout << "#pos_i pos_j original_score apc_score lrs_score\n";
  for (int i = 0; i < M.rows(); i++) {
    for (int j = i + 1; j < M.cols(); j++) {
      fout << i + 1 << " " << j + 1 << " " << M(i, j) << " " << apc(i, j) << " "
           << S(i, j) << std::endl;
    }
  }

  return 0;
}

}  // namespace colors
