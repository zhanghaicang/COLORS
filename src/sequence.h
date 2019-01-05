#ifndef COLORS_SEQUENCE_H_
#define COLORS_SEQUENCE_H_

#include <cstdio>
#include <string>
#include <vector>

#include "type.h"

namespace colors {

int read_msa(const char* file_path, std::vector<std::string>& msa);
unsigned char aatoi(const unsigned char aa);

int calc_seq_weight(const std::vector<std::string>& msa, const double seq_id,
                    Vector& weight);
double calc_seq_sim(const std::string& seq1, const std::string& seq2);
}  // namespace colors

#endif
