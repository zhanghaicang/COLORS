#include "score.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace colors {

int calc_matrix(const std::vector<std::string>& msa, const Vector& weight,
                Matrix& mi, Matrix& omes, Matrix& cov) {
  int nrow = msa.size();
  int ncol = msa[0].size();

  const int pseudo_count = 1.0;

  Matrix pa(ALPHA, ncol);
  Matrix pab(ALPHA, ALPHA);

  double neff = weight.sum();
  pa.setConstant(pseudo_count);
  for (int i = 0; i < ncol; i++) {
    auto col = pa.col(i);
    for (int j = 0; j < nrow; j++) {
      col(msa[j][i]) += weight[j];
    }
    col /= 1.0 * pseudo_count * ALPHA + neff;
  }

  for (int i = 0; i < ncol; i++) {
    for (int j = i; j < ncol; j++) {
      if (i == j) {
        pab.setZero();
        pab.diagonal() = pa.col(i);
      } else {
        pab.setConstant(pseudo_count * 1.0 / ALPHA);
        for (int k = 0; k < nrow; k++) {
          pab(msa[k][i], msa[k][j]) += weight(k);
        }
        pab /= 1.0 * pseudo_count * ALPHA + neff;
        mi(i, j) = mi(j, i) = calc_mi(pab, pa.col(i), pa.col(j));
        omes(i, j) = omes(j, i) = calc_omes(pab, pa.col(i), pa.col(j));
        cov(i, j) = cov(j, i) = calc_cov(pab, pa.col(i), pa.col(j));
      }
    }
  }

  return 0;
}

double calc_mi(const Matrix& pab, const Vector& pa, const Vector& pb) {
  return (pab.array() *
          Eigen::log((pab.array() / (pa * pb.transpose()).array())))
      .sum();
}

double calc_omes(const Matrix& pab, const Vector& pa, const Vector& pb) {
  return (pab.array() - (pa * pb.transpose()).array()).square().sum();
}

double calc_cov(const Matrix& pab, const Vector& pa, const Vector& pb) {
  return (pab.array() - (pa * pb.transpose()).array()).abs().sum();
}

int calc_apc(const Matrix& m, Matrix& apc) {
  auto& row_mean = (m.rowwise().sum() - m.diagonal()) / (m.cols() - 1);
  apc = (m - row_mean * row_mean.transpose()) / row_mean.mean();

  return 0;
}

}  // namespace colors
