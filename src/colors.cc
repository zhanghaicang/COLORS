#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "lrs.h"
#include "score.h"
#include "sequence.h"
#include "type.h"

using namespace colors;

int main(int argc, const char** argv) {
  if (argc != 3) {
    std::cout << "Usage: <1>msa-file <2>output-prefix\n";
    return -1;
  }

  const char* msa_path = argv[1];
  const char* output_prefix = argv[2];
  std::vector<std::string> msa;
  read_msa(msa_path, msa);
  int nrow = msa.size(), ncol = msa[0].size();

  Vector weight(nrow);
  calc_seq_weight(msa, 0.8, weight);

  std::cout << "nrow= " << nrow << " ncol= " << ncol
            << " neff= " << weight.sum() << std::endl;

  Matrix mi(ncol, ncol);
  Matrix omes(ncol, ncol);
  Matrix cov(ncol, ncol);
  calc_matrix(msa, weight, mi, omes, cov);

  std::string mi_output = std::string(output_prefix) + ".mi.lrs";
  calc_lrs(0.13, mi, mi_output.c_str());

  std::string omes_output = std::string(output_prefix) + ".omes.lrs";
  calc_lrs(0.26, omes, omes_output.c_str());

  std::string cov_output = std::string(output_prefix) + ".cov.lrs";
  calc_lrs(0.1, cov, cov_output.c_str());

  return 0;
}
