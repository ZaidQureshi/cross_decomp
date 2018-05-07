#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <cstdint>
#include <random>
#include <algorithm>

using namespace std;

#define NUM_NODES 500

void initMat(bool*& adjMat) {
  adjMat = new bool[NUM_NODES*NUM_NODES];
  std::fill(adjMat, adjMat + NUM_NODES * NUM_NODES ,0); // Ketan Date: Initialized adjMat to 0.
  //#pragma omp parallel for
  for(int i =0 ; i < NUM_NODES; i++)
	  adjMat[i * NUM_NODES + i] = true; // Ketan Date: Added self loop
  
}

void readFile(bool*& adjMat, const char* filename) {
  uint64_t u, v, w;
  ifstream initFile(filename);
  while (initFile >> u >> v >> w) {
	  adjMat[(u-1) * NUM_NODES + (v-1)] = true;  // Ketan Date: Made explicit to be sure.
	  adjMat[(v-1) * NUM_NODES + (u-1)] = true;
//    if (!(adjMat[(v-1) * NUM_NODES + (u-1)]))
      
  }

  
}

void freeMat(bool*& adjMat) {
  delete adjMat;
}

void printOrigMat(bool* adjMat) {

	for (size_t i = 0; i < (NUM_NODES*NUM_NODES); i++) {
		if (i != 0 && ((i % (NUM_NODES)) == 0))
		cout << endl;
		cout << " " << (int) adjMat[i];
	}
}
// void printMatPointsWithPartitions(bool*& adjMat,
// 				  vector<vector<uint64_t>>&row_parts,
// 				  vector<vector<uint64_t>>&col_parts,
// 				  const uint64_t& num_rows,
// 				  const uint64_t& num_cols,
// 				  const uint64_t& v) {
//   uint64_t i, j;
//   for (j = 0; j < v; j++) {
//     cout <<
//   }
// }

void printParts(vector<vector<uint64_t>>& row_parts, vector<vector<uint64_t>>& col_parts) {
  
  for (size_t i = 0; i < row_parts.size(); i++) {
    cout << "ROW PARTITION\t" << i << endl;
    vector<uint64_t>& curPart = row_parts[i];
    //sort(curPart.begin(), curPart.end());
    for (size_t j = 0; j < curPart.size(); j++)
      cout << "\t" << curPart[j] << endl;
  }
  for (size_t i = 0; i < col_parts.size(); i++) {
    cerr << "COL PARTITION\t" << i << endl;
    vector<uint64_t>& curPart = col_parts[i];
    //sort(curPart.begin(), curPart.end());
    for (size_t j = 0; j < curPart.size(); j++)
      cerr << "\t" << curPart[j] << endl;
  }
}

void initParts(vector<vector<uint64_t>>& row_parts,
	       vector<vector<uint64_t>>& col_parts,
	       const uint64_t& v,
	       const uint64_t& num_rows,
	       const uint64_t& num_cols) {
  random_device rd;
  mt19937 mt(rd());
  uniform_int_distribution<int> dist(0, v-1);
  //#pragma omp parallel for
  for (size_t i = 0; i < num_rows; i++) {
    uint64_t part = dist(mt);
    row_parts[part].push_back(i);
    col_parts[part].push_back(i);
  }
  //for (size_t i = 0; i < num_cols; i++) {
    //uint64_t part = dist(mt);
    //row_parts[part].push_back(i);
    //col_parts[0].push_back(i);
    //}
}

void initWeights(float*& alpha, float*& beta,
		 const uint64_t& num_rows, const uint64_t& num_cols) {
	alpha = new float[num_rows];
	beta = new float[num_cols];
	#pragma omp parallel for
	for (size_t i = 0; i < num_rows; i++)
		alpha[i] = 1;
	#pragma omp parallel for
	for (size_t i = 0; i < num_cols; i++)
		beta[i] = 1;
}

void repartition(bool*& adjMat, float*& alpha, float*& beta,
		 vector<vector<uint64_t>>& row_parts, vector<vector<uint64_t>>& col_parts,
		 const uint64_t& num_rows, const uint64_t& num_cols, const uint64_t& v,
		 const float h, const float lambda) {

  //build new partition on columns
  #pragma omp parallel for
  for (size_t j = 0; j < num_cols; j++) {
    vector<float> y_j(v);
    //for each partition
    
    for (size_t r = 0; r < v; r++) {
      float y_j_r = beta[j];
      //first part of objective function
      float temp1 = h;
      vector<uint64_t>& curRowPart = row_parts[r];
      size_t curRowPartSize = curRowPart.size();
      float sum  = 0;
      //for each row in the partition
      for (size_t i = 0; i < curRowPartSize; i++) {
	uint64_t row = curRowPart[i];
	sum += alpha[row] * pow((float)adjMat[(row * num_cols) + j ], (float)lambda);
      }
      temp1 *= sum;
      //second part of objective function
      float temp2 = (1 - h);
      sum = 0;
      for (size_t i = 0; i < num_rows; i++) {
	if (find(curRowPart.begin(), curRowPart.end(), i) == curRowPart.end()) {
	  sum += alpha[i] * pow(1 - (float) adjMat[i * num_cols + j], (float)(1/lambda));
	}
      }
      temp2 *= sum;
      y_j_r *= (temp1 + temp2);
      y_j[r] = y_j_r;
    }
    //find partition with max objective function
    size_t maxIdx = 0;
    for (size_t x = 1; x < v; x++) {
      if (y_j[maxIdx] < y_j[x])
	maxIdx = x;
    }
    //remove column from previous partition
    //#pragma omp parallel for
    for (size_t x = 0; x < v; x++) {
      vector<uint64_t>::iterator it1 = col_parts[x].end();
      vector<uint64_t>::iterator it2 = remove(col_parts[x].begin(), col_parts[x].end(), j);
      if (it1 != it2)
	col_parts[x].pop_back();
	
    }
    //assign to new partition
    col_parts[maxIdx].push_back(j);
  }
  //build new partition on rows
  #pragma omp parallel for
  for (size_t i = 0; i < num_rows; i++) {
    vector<float> y_i(v);
    //for each partition
    //#pragma omp parallel for
    for (size_t r = 0; r < v; r++) {
      float y_i_r = alpha[i];
      //first part of the objective function
      float temp1 = h;
      vector<uint64_t>& curColPart = col_parts[r];
      size_t curColPartSize = curColPart.size();
      float sum = 0;
      //for each column in the partition
      for (size_t j = 0; j < curColPartSize; j++) {
	uint64_t col = curColPart[j];
	sum += beta[col] * pow(adjMat[(i * num_cols) + col], lambda);
	
      }
      temp1 *= sum;
      //second part of the objective function
      float temp2 = (1 - h);
      sum = 0;
      for (size_t j = 0; j < num_cols; j++) {
	if (find(curColPart.begin(), curColPart.end(), j) == curColPart.end()) {
	  sum += beta[j]  * pow(1 - adjMat[i * num_cols + j], (1/lambda));
	}
      }
      temp2 *= sum;
      y_i_r *= (temp1 + temp2);
      y_i[r] = y_i_r; // Ketan Date: Corrected from y_i[v] = y_i_r;
      
    }
    //find partition with max objective function
    size_t maxIdx = 0;
    for (size_t x = 1; x < v; x++) {
      if (y_i[maxIdx] < y_i[x])
	maxIdx = x;
    }
    //remove row from previous partition
    //#pragma omp parallel for
    for (size_t x = 0; x < v; x++) {
      vector<uint64_t>::iterator it1 = row_parts[x].end();
      vector<uint64_t>::iterator it2 = remove(row_parts[x].begin(), row_parts[x].end(), i);
      if (it1 != it2)
	row_parts[x].pop_back();
    }
    //assign to new partition
    row_parts[maxIdx].push_back(i);
  }
}

int main(int argc, char** argv) {
	bool* adjMat;
	float h;
	float lambda;
	uint64_t v;
	float* alpha;
	float* beta;
	uint64_t num_rows;
	uint64_t num_cols;
	//parse arguments
	if (argc < 5)
	    exit(1);
	else {
	    h = atof(argv[2]);
	    lambda = atof(argv[3]);
	//    v = _strtoull(argv[4], nullptr, 10);
		v = 10; // Ketan Date: Fixed value of v for testing. strtoull not recognized in VS 2012
	}
	vector<vector<uint64_t>> row_parts(v);
	vector<vector<uint64_t>> col_parts(v);
	num_rows = NUM_NODES;
	num_cols = NUM_NODES;
	//init adjacency matrix
	initMat(adjMat);
	//init partitions
	initParts(row_parts, col_parts, v, num_rows, num_cols);
	//init weights
	initWeights(alpha, beta, num_rows, num_cols);
	
	//read init graph
	readFile(adjMat, argv[1]);
	  
	//print orig mat
//	  printOrigMat(adjMat);
	//repartition
	//printParts(row_parts, col_parts);
	for (uint64_t i = 0; i < 10; i++) {
	    repartition(adjMat, alpha, beta, row_parts, col_parts,
			num_rows, num_cols, v,
			h, lambda);
	    
	    //cout << "********************************" << endl;
	    //	printParts(row_parts, col_parts);
	}
	//check convergence

	//print matrix points with partitions
	// printMatPointsWithPartitions(adjMat, row_parts, col_parts,
	// 			       num_rows, num_cols, v);
	  
	//free adjacency matrix
	freeMat(adjMat);
	  
}
