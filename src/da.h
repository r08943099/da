#ifndef DA_H
#define DA_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "node.h"
using namespace std;
#define citysize 14
#define iteration 50000
struct DA
{
    Node                _nodeArray[citysize] = {};
    string              _problemName;
    double              _distance_matrix[citysize][citysize] = {};
    bool                _qubit_matrix[citysize][citysize] = {};
    double              _A = 1, _B = 2 , _C = 2; //BC need add
    double              _beta = 0,_min_beta = 0.01, _max_beta = 2;
    double              _E_off = 0;
    double              _E_off_increment = 10;     
};

void parseInput2da(fstream& inFile, struct DA*);
void calculate_distnace(struct DA*);
bool calculate_delta(struct DA*);
// double calculate_energy(struct DA* da);
// double calculate_delta_energy(struct DA* da, int i, int j)
bool replica(struct DA*);

#endif // DA_H

