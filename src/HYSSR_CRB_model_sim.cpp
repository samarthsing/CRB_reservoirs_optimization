#include <bits/stdc++.h>
#include <numeric>
#include <mpi.h>
#include "HYSSR_CRB_model.cpp"

using namespace std;



int main(int argc, char* argv[])
{
         int world_rank; // the rank of the process
         int world_size; // number of processes
         
         vector<vector<double> > vars_matrix=utils::readMatrixFromDataFile("policy_to_be_simulated.txt");
         cout<<vars_matrix.size()<<endl;
         
         
         MPI_Init(NULL, NULL);      // initialize MPI environment
         
         MPI_Comm_size(MPI_COMM_WORLD, &world_size);

         
         MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
         printf("I evaluated rank = %d\n",world_rank);
         
         
         //vector<double> vars=utils::loadVectorFromDataFile("min_spill_pareto.txt");
          
         //vector<vector<double> > objs_matrix=vector<vector<double> >(vars_matrix.size());
         vector<double> objs(4,0);
         vector<double> constrs(0,0);
          
         vector<double> vars;
         vars = vars_matrix[world_rank];
           
         MyClass obj1=MyClass();
         //Fcrps obj2=Fcrps(vars);
         
         double* var_array = &vars[0];
         double* obj_array = &objs[0];
         double* constrs_array = &constrs[0];
         obj1.evaluate(var_array,obj_array, constrs_array,world_rank,1);
         
         
        MPI_Finalize();
         
         return 0;

}