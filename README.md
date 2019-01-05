## Overview
COLORS (improving COntact prediction using LOw-Rank and Sparse matrix decomposition) is an approach to seperate the signal from the background in the correlation analysis for protein contacts prediction.  In our approach, a correlation matrix was decomposed into two components, i.e., a low-rank component representing background correlations, and a sparse component representing true correlations. Finally the residue contacts were inferred from the sparse component of correlation matrix. 


## Citation
Haicang Zhang, Yujuan Gao, Minghua Deng, Chao Wang, Jianwei Zhu, Shuai Cheng Li, Wei-Mou Zheng, Dongbo Bu. Improving residue–residue contact prediction via low-rank and sparse decomposition of residue correlation matrix. Biochemical and biophysical research communications 472.1 (2016): 217-222.

## Build
1. Install the Eigen library. Plese refere to https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html for details.
2. cd COLORS/src; make

## Usage
./colors input-MSA ouptut-prefix

