#ifndef DA_H
#define DA_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "node.h"

using namespace std;
#define citysize 14
#define iteration 5   //for each setting
#define repeat    10  //for repeat this code
struct DA
{
    Node                _nodeArray[citysize] = {};
    string              _problemName;
    double              _distance_matrix[citysize][citysize] = {};
    bool                _qubit_matrix[citysize][citysize] = {};
    double              _A = 1, _B, _C; //BC need add
    double              _E_off = 0;
    double              _E_off_increment = 2;     
    double              _r = 0.999;
    double              _beta;
    double              _best_energy;
    bool                _best_qubit_matrix[citysize][citysize] = {};
};

void parseInput2da(fstream& inFile, struct DA*);
void reset(struct DA* da);
void calculate_distnace(struct DA*);
double calculate_energy(struct DA*);
bool calculate_delta(struct DA*);
// double calculate_energy(struct DA* da);
// double calculate_delta_energy(struct DA* da, int i, int j)
bool replica(struct DA*);

#endif // DA_H

