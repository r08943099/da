#include "da.h"
#include <time.h>
#include "math.h"
//#include "gnuplot-iostream.h"

int main(int argc, char** argv)
{
    srand( time(NULL) );
    fstream input, output;
    if (argc == 3) {
        input.open(argv[1], ios::in);
        output.open(argv[2], ios::out);
        if (!input) {
            cerr << "Cannot open the input file \"" << argv[1]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./da <input file> <output file>" << endl;
        exit(1);
    }
    DA da;
    parseInput2da(input, &da);
    // da._B = 2000/da._beta;
    // da._C = 2000/da._beta;
    //da._beta *= da._r;
    //double rep_best_energy_list[repeat] = {};
    double sum_best_energy = 0;
    double rep_best_energy = 10000;
    bool   rep_best_qubit_matrix[citysize][citysize] = {};
    for(int rep = 0; rep < repeat; rep++){//repeat this code for x times    
        cout << "=============================== Repeat this code : " << rep << " ================================" <<endl;
        reset(&da);
        while(da._beta <= 2){
            da._beta /= da._r;
            for(int i = 0; i < iteration; i++){
                //cout << "=============================== The trial B: " << i << " ================================" <<endl;
                if(!replica(&da)){
                    cout << "Cannot escape local Minimun" << endl;
                    break;
                }
            }
        }
        //rep_best_energy_list[repeat] = da._best_energy;
        sum_best_energy += da._best_energy;
        if(da._best_energy < rep_best_energy){
            rep_best_energy = da._best_energy;
            for(int i = 0; i < citysize; i++){
                for(int j = 0; j < citysize; j++){
                    rep_best_qubit_matrix[i][j] = da._best_qubit_matrix[i][j];
                }       
            }
        } 
    }      
  
    //print the benchmark information
    cout << "==== Print Result ====" << endl;
    cout << "-------------------------" << endl;
    cout << "The Best Cost : " << rep_best_energy << endl;
    cout << "The Average Cost : " << sum_best_energy/repeat << endl;
    cout << "The Best Qubit : " << endl;
    for(int i = 0; i < citysize; i++){
        for(int j = 0; j < citysize; j++){
            cout << rep_best_qubit_matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << "-------------------------" << endl;
    cout << "=========================================" << endl;
    
    //output file
    for(int i = 0; i < citysize; i++){
        for(int j = 0; j < citysize; j++){
            if(rep_best_qubit_matrix[i][j] == 1){
                output << da._nodeArray[j]._x << " " << da._nodeArray[j]._y << endl;
            }
            
        }
    }

    output.close();
    
    return 0;
}