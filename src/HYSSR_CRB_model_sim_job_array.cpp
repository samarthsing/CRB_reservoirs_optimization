#include <bits/stdc++.h>
#include <numeric>
#include <mpi.h>
#include "HYSSR_CRB_model.cpp"

using namespace std;



int main(int argc, char* argv[])
{
         
         vector<vector<double> > vars_matrix=utils::readMatrixFromDataFile("policies_to_be_simulated_full.txt");
         cout<<vars_matrix.size()<<endl;
         int world_rank=atoi(argv[1]);
         printf("I evaluated rank = %d\n",world_rank);
         vector<double> objs(6,0);
         vector<double> constrs(1,0);
         vector<double> vars;
         vars = vars_matrix[world_rank];
         cout<<"vars.size() "<<vars.size()<<endl;
         MyClass obj1=MyClass();
         double* var_array = &vars[0];
         double* obj_array = &objs[0];
         double* constrs_array = &constrs[0];
         obj1.evaluate(var_array,obj_array, constrs_array,world_rank,1);
         
         return 0;

}