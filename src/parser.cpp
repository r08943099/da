// /*========== Parse a Benchmarks file ==========*/
// /* The format is:
//     NAME:<problemName>
//     DIMENSION:<nodeNum>
//     NODE_COORD_SECTION:
//         <nodeIdx> <x> <y>
// */
// #include "da.h"
// string& ClearAllSpace(string &str)
// {
//     int index = 0;
//     if( !str.empty())
//     {
//         while( (index = str.find(' ',index)) != string::npos)
//         {
//             str.erase(index,1);
//         }
//     }

//     return str;
// }
// void parseInput2da(fstream& inFile, struct DA* da)
// {
//     string str;
//     while(getline(inFile,str)){
//         ClearAllSpace(str);
//         int idx = str.find(":");
//         string dataType = str.substr(0,idx);
//         string data = str.substr(idx+1); 
//         if(dataType == "NAME") da -> _problemName = data;
//         if(dataType == "NODE_COORD_SECTION")break;
//     }
//     for(int i = 0; i < citysize; i++){
//         inFile >> str;//neglect name
//         inFile >> str;
//         da -> _nodeArray[i]._x = stod(str);
//         inFile >> str;
//         da -> _nodeArray[i]._y = stod(str);
//     }
//     //initial patameter
//     for(int i = 0; i < citysize; i++){
//         for(int j=0; j < citysize; j++){
//             if(i == j)da -> _qubit_matrix[i][j] = 0;
//         }
//     }
//     // da -> _qubit_matrix[0][0] = 1;da -> _qubit_matrix[0][1] = 1;da -> _qubit_matrix[0][2] = 1;
//     // da -> _qubit_matrix[1][0] = 1;da -> _qubit_matrix[1][1] = 0;da -> _qubit_matrix[1][2] = 0;
//     // da -> _qubit_matrix[2][0] = 0;da -> _qubit_matrix[2][1] = 0;da -> _qubit_matrix[2][2] = 0;
//     // double longest_distnace = 0;
//     // calculate_distnace(da);
//     // for(int i = 0; i < citysize; i++){
//     //     for(int j=0; j < citysize; j++){
//     //         if(longest_distnace < da -> _distance_matrix[i][j]) longest_distnace = da -> _distance_matrix[i][j];
//     //     }
//     // }
//     // da -> _B = longest_distnace*2;
//     // da -> _C = longest_distnace*2;
//     // cout << " The original energy  = " << calculate_energy(da) << endl;
//     // //print the benchmark information
//     // cout << "==== Print the benchmark information ====" << endl;
//     // cout << "NAME : " << da -> _problemName << endl;
//     // cout << "DIMENSION : " << citysize << endl;
//     // cout << "NODE_COORD_SECTION" << endl;
//     // cout << "-------------------------" << endl;
//     // for(int i = 0; i < citysize; i++){
//     //     cout << "NodeIdx : " << i << " ";
//     //     cout << "X : "       << da -> _nodeArray[i]._x       << " ";
//     //     cout << "Y : "       << da -> _nodeArray[i]._y       << endl;
//     // }
//     // // cout <<  "Distance matrix: "<< endl;
//     // // for(int i = 0; i < citysize; i++){
//     // //     for(int j=0; j < citysize; j++){
//     // //         cout << da -> _distance_matrix[i][j] << " ";
//     // //     }
//     // //     cout << endl;
//     // // }
//     // cout << "-------------------------" << endl;
//     // cout << "=========================================" << endl;
//     return;
// }