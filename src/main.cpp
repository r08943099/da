#include "da.h"
#include <time.h>

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
    for(int i = 0; i < iteration; i++){
        cout << "The trial B: " << i << endl;
        da._beta += (da._max_beta - da._min_beta) / iteration; //-> pointer use
        if(!replica(&da)){
            cout << "Cannot escape local Minimun" << endl;
            break;
        }
    }
    // //print the benchmark information
    // cout << "==== Print Result ====" << endl;
    // cout << "-------------------------" << endl;
    // cout << "Qubit : " << endl;
    // for(int i = 0; i < citysize; i++){
    //     for(int j = 0; j < citysize; j++){
    //         cout << da._qubit_matrix[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "-------------------------" << endl;
    // cout << "=========================================" << endl;
    return 0;
}