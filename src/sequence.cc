#include "sequence.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace colors {

int read_msa(const char* file_path, std::vector<std::string>& msa) {
  std::ifstream fin;
  fin.open(file_path);
  if (fin.fail()) {
    std::cerr << "Failed in open file " << file_path << std::endl;
    return -1;
  }

  std::string line;
  while (getline(fin, line, '\n')) {
    if (line.size() == 0 || line[0] == '>') {
      continue;
    }
    std::for_each(line.begin(), line.end(), [](char& c) { c = aatoi(c); });
    msa.push_back(line);
  }
  fin.close();

  return 0;
}

int calc_seq_weight(const std::vector<std::string>& msa, const double seq_id,
                    Vector& weight) {
  int nrow = msa.size();
  weight.setOnes();

  for (int i = 0; i < nrow; i++) {
    for (int j = i + 1; j < nrow; j++) {
      double sim = calc_seq_sim(msa[i], msa[j]);
      if (sim >= seq_id) {
        weight[i] += 1.0;
        weight[j] += 1.0;
      }
    }
  }
  weight = weight.cwiseInverse();

  return 0;
}

double calc_seq_sim(const std::string& seq1, const std::string& seq2) {
  double sim = 0.0;
  for (int i = 0; i < seq1.size(); i++) {
    sim += seq1[i] == seq2[i];
  }

  return sim / seq1.size();
}

unsigned char aatoi(const unsigned char aa) {
  unsigned char id;
  switch (aa) {
    case '-':
      id = 0;
      break;
    case 'A':
      id = 1;
      break;
    case 'C':
      id = 2;
      break;
    case 'D':
      id = 3;
      break;
    case 'E':
      id = 4;
      break;
    case 'F':
      id = 5;
      break;
    case 'G':
      id = 6;
      break;
    case 'H':
      id = 7;
      break;
    case 'I':
      id = 8;
      break;
    case 'K':
      id = 9;
      break;
    case 'L':
      id = 10;
      break;
    case 'M':
      id = 11;
      break;
    case 'N':
      id = 12;
      break;
    case 'P':
      id = 13;
      break;
    case 'Q':
      id = 14;
      break;
    case 'R':
      id = 15;
      break;
    case 'S':
      id = 16;
      break;
    case 'T':
      id = 17;
      break;
    case 'V':
      id = 18;
      break;
    case 'W':
      id = 19;
      break;
    case 'Y':
      id = 20;
      break;
    default:
      id = 0;
  }
  return id;
}

}  // namespace colors
