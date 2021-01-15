#include "da.h"
#include "math.h"

void reset(struct DA* da){
    for(int i = 0; i < citysize; i++){
        for(int j=0; j < citysize; j++){
            da ->_qubit_matrix[i][j] = 1;
        }
    }
    da -> _beta = 0.02;   
    da -> _best_energy = 10000;     
}

void calculate_distnace(struct DA* da)
{
    for(int i = 0; i < citysize; i++){
        for(int j = 0; j < citysize; j++){
            da -> _distance_matrix[i][j] = 
            sqrt(pow((da -> _nodeArray[i]._x - da -> _nodeArray[j]._x),2) + pow((da -> _nodeArray[i]._y - da -> _nodeArray[j]._y),2));        
        }
    }
}

int penalty_funciotn(int x){
    x = (x == 0)?10:x-1;
    //x = x-1;
    return pow(x,2);
}

double calculate_energy(struct DA* da)
{   
    //i order
    //k city
    double H = 0, HA = 0, HB = 0, HC=0;
    //HA,HC
    for(int k = 0; k < citysize; k++){    
        int Bi = 0, Ci = 0;    
        for(int i = 0; i < citysize; i++){
            if(da -> _qubit_matrix[i][k]){
                Ci ++;
                if(i < citysize - 1){
                    for(int l = 0; l < citysize; l++){
                        HA += da -> _distance_matrix[k][l] * da -> _qubit_matrix[i+1][l];
                    }
                }
                else{ //last city to first city
                    for(int l = 0; l < citysize; l++){
                        HA += da -> _distance_matrix[k][l] * da -> _qubit_matrix[0][l];
                    }
                }

            }
            if(da -> _qubit_matrix[k][i]) Bi++;
        }
        HB += penalty_funciotn(Bi);
        HC += penalty_funciotn(Ci);
    }   
    H = da -> _A*HA + da -> _B*HB + da -> _C*HC;
    return H;
}

double calculate_delta_energy(struct DA* da, int i, int j)
{
    double H = 0, HA = 0 , HB = 0, HC = 0;
    int sgn;
    int Bi = 0, Ci = 0; 
    sgn = (da ->_qubit_matrix[i][j])?-1:1; //sign
    for(int l = 0; l < citysize; l++){ 
        for(int m = i - 1; m < i + 2; m+=2){  
            if(m >= 0 && m < citysize)  HA += da -> _distance_matrix[j][l] * da -> _qubit_matrix[m][l];
        }
        if(da -> _qubit_matrix[i][l]){
            Bi ++;
        }
        if(da -> _qubit_matrix[l][j]){
            Ci ++;
        }
    }
    HB = penalty_funciotn(Bi+sgn) - penalty_funciotn(Bi); 
    HC = penalty_funciotn(Ci+sgn) - penalty_funciotn(Ci);   
    H = sgn * da -> _A*HA + da -> _B*HB + da -> _C*HC;
    return H;
}


bool ADB(struct DA* da, double delta_energy)
{
    delta_energy = delta_energy - da -> _E_off; // add or minus ?
    //an ramdom number[0,1)
    double random0to1 = (double) rand() / (RAND_MAX + 1.0);
    double probabilty = exp(-delta_energy/ da -> _beta);
    double acceptance_probability;  
    if(probabilty < 1) acceptance_probability = probabilty;
    else acceptance_probability = 1;
    //cout << " delta_energy: " << delta_energy << " random0to1:" << random0to1 << " probabilty:" << probabilty <<  " acceptance_probability: " << acceptance_probability << endl;
    if(acceptance_probability >= random0to1) return 1; //mean accept
    else return 0;
    // return delta_energy < 0;
}
bool random_choose_flip(bool candicate[citysize][citysize], int* flipx, int* flipy){
    int candicatecount = 0;
    for(int i = 0; i < citysize; i++){
        for(int j = 0; j < citysize; j++){
            if(candicate[i][j]) candicatecount++;  
        }
    }
    if(candicatecount == 0) return false; //no cadicate;        
    int ramdom_candicate = rand() % candicatecount;
    //cout << "ramdom_candicate" << ramdom_candicate << endl;
    candicatecount = 0;
    for(int i = 0; i < citysize; i++){
        for(int j = 0; j < citysize; j++){
            if(candicate[i][j]){  
                if(candicatecount++ == ramdom_candicate){
                    *flipx = i;
                    *flipy = j;
                    return true;  
                }
            }    
        }
    }
    return false;     
}  

void find_best_energy(struct DA* da, double energy){
    if(energy < da -> _best_energy){
        da -> _best_energy = energy;
        for(int i = 0; i < citysize; i++){
            for(int j = 0; j < citysize; j++){
                 da -> _best_qubit_matrix[i][j] = da -> _qubit_matrix[i][j]; 
            }
        }
    }
}

//parallel
bool calculate_delta(struct DA* da)
{
    bool flipResult = 0;
    int flipx, flipy;
    da -> _E_off = 0;
    while(!flipResult){
        bool candicate[citysize][citysize] = {};
        //cout << "Candicate: " << endl; 
        for(int i = 0; i < citysize; i++){
            for(int j = 0; j < citysize; j++){
                double delta_energy = calculate_delta_energy(da, i, j);//many candicate
                //cout << i<< ":delta energy" << delta_energy << endl; 
                //cout << ADB(da, delta_energy) << endl; 
                //cout<< "Location : "<< "i = " << i << " j=" << j << " is : " << endl;
                candicate[i][j] = ADB(da, delta_energy);  
                //cout << "The ADB result : "<<candicate[i][j] << endl;//if both less and more both enter 
            }
            //cout << endl;
        }
        //choose one to flip;
        flipResult = random_choose_flip(candicate, &flipx, &flipy);
        da -> _E_off += da -> _E_off_increment; 
        if(da ->_E_off > 100000) {
            //cout << "Eoff = "<< da ->_E_off << endl;
            return false;
        }
    }
    //cout << "The beta :" << da -> _beta << endl;   
    da ->  _qubit_matrix[flipx][flipy] = !da ->  _qubit_matrix[flipx][flipy];
    //print the benchmark information
    //cout << "==== Print Result ====" << endl;
    //cout << "-------------------------" << endl;
    double energy = calculate_energy(da); //del after
    find_best_energy(da, energy); //might edit
    //cout << "After flip Qubit : " << endl;
    //cout << "Total Energy" << energy << endl;
    //cout << "The Best Energy" << da -> _best_energy << endl;
    //cout << "The Current Qubit :" << endl;
    // for(int i = 0; i < citysize; i++){
    //     for(int j = 0; j < citysize; j++){
    //         cout << da -> _qubit_matrix[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "The Best Qubit :" << endl;
    // for(int i = 0; i < citysize; i++){
    //     for(int j = 0; j < citysize; j++){
    //         cout << da -> _best_qubit_matrix[i][j] << " ";
    //     }
    //     cout << endl;
    // }    
    //cout << "-------------------------" << endl;
    //cout << "=========================================" << endl;  
    return true;     
}

bool replica(struct DA* da){
    return calculate_delta(da); 
}