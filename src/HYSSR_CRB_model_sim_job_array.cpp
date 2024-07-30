#include <bits/stdc++.h>
#include <numeric>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
#include "HYSSR_CRB_model.cpp"

using namespace std;



int main(int argc, char* argv[])
{
         
         vector<vector<double> > vars_matrix=utils::readMatrixFromDataFile("policies_to_be_simulated_full.txt");
         cout<<vars_matrix.size()<<endl;
         cout<<"argv[1] "<<argv[1]<<endl;
         int world_rank=atoi(argv[1]);
         //int master_rank=0;
         //int master_rank=atoi(argv[1])%9;
         int master_rank=2;
         //int world_rank=0;
         //int master_rank=2;
         printf("I evaluated rank = %d\n",world_rank);
         vector<double> objs(4,0);
         vector<double> constrs(1,0);
         vector<double> vars;
         vars = vars_matrix[world_rank];
         cout<<"vars.size() "<<vars.size()<<endl;
         MyClass obj1=MyClass();
         double* var_array = &vars[0];
         double* obj_array = &objs[0];
         double* constrs_array = &constrs[0];
         cout<<"world_rank "<<world_rank<<endl;
         cout<<"master_rank "<<master_rank<<endl;
         obj1.evaluate(var_array,obj_array, constrs_array,world_rank,master_rank);
         
         return 0;

}