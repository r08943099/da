#ifndef DA_H
#define DA_H

#include <iostream>
#include <fstream>
#include "node.h"
// #include "ap_fixed.h"

using namespace std;
#define citysize 4
#define iteration 100  //for each setting
#define repeat    10  //for repeat this code
#define replicaNum 1000
#define beta 0.00001
#define cooling_rate 0.85
#define _E_off_increment 1

// struct Replica
// {
//     double _beta;
//     bool   _qubit_matrix[citysize][citysize] = {};
// };

struct Replica
{
    double _beta;
    bool   _qubit_matrix[citysize][citysize];
    //bool   _best_qubit_matrix[citysize][citysize]; 
    double _energy;
};

struct DA
{   
    Node                _nodeArray[citysize];
    double              _distance_matrix[citysize][citysize];
    bool                _current_qubit_matrix[citysize][citysize];
    double              _A;
    double              _E_off;
    //double              _E_off_increment;     
    double              _beta;
    double              _best_energy;
    bool                _best_qubit_matrix[citysize][citysize];
    Replica             _replicaArray[replicaNum];
};

void DigitalAnnealer(struct DA*);
// void reset(struct DA* da);
// void calculate_distnace(struct DA*);
// double calculate_energy(struct DA*);
// bool calculate_delta(struct DA*);
// // double calculate_energy(struct DA* da);
// // double calculate_delta_energy(struct DA* da, int i, int j)
// bool replica(struct DA*);

#endif // DA_H

