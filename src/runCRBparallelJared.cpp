//============================================================================
// Name        : mainParallel_lakecomo.cpp
// Author      : MatteoG
// Version     :
// Copyright   : Your copyright notice
//============================================================================

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <fstream>
#include "HYSSR_CRB_model.cpp"
#include "../borg_jared/borgmm.h"

using namespace std;

#define PI 3.14159265

void model_wrapper(double *vars, double *objs, double *constr, int *evals, int *isle, int *worker)
{
    // run the model. Change monte carlo evaluations here.
    //cout<<"nfe inside model_wrapper is "<<*evals<<endl;
    unsigned int N = *evals;
    //cout<<"unsigned nfe inside model_wrapper is "<<N<<endl;
    unsigned int mstr=*isle;
    MyClass a=MyClass();
    a.evaluate(vars, objs, constr,N,mstr);
}


int main(int argc, char* argv[]) {

    //feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);

    // Red River Problem
    //int nvars = 281;
    //int nvars=26;
    int nvars = 144;
    int nobjs = 4;
    int nconstrs = 1;
    //BORG_Debug_on();
    // setting random seed
    unsigned int seed = atoi(argv[1]);
    srand(seed);
    //int NFE = atoi(argv[2]);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Changes for MPI version start here
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    BORG_Algorithm_ms_startup(&argc, &argv);
    BORG_Algorithm_ms_max_time(stod(argv[2]));
    //BORG_Algorithm_ms_max_evaluations(atoi(argv[2])); // Choose total NFE here
    BORG_Algorithm_ms_islands(atoi(argv[3]));

    //Enable global Latin hypercube initialization to ensure each island
    // gets a well sampled distribution of solutions.
    BORG_Algorithm_ms_initialization(INITIALIZATION_LATIN_GLOBAL);

    
    
    BORG_Algorithm_output_frequency(200); // print archive every NFE/100 evaluations

    // Define the problem with decisions, objectives, constraints, and the evaluation function
    BORG_Problem problem = BORG_Problem_create(nvars, nobjs, nconstrs, model_wrapper);

    // Set all the parameter bounds and epsilons
    unsigned int j=0;
    
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++;// constant for non-convex RBF
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++;// constant for non-convex RBF
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++;// constant for non-convex RBF
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++;// constant for non-convex RBF
    for(unsigned int i=0; i<8; i++){ // loop through 15 RBFs
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 1
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 1
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 2
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 2
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 3
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 3
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 4
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 4
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 5
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 5
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for sine term
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for sine term
        BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for cosine term
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for cosine term
        //BORG_Problem_set_bounds(problem, j, -1.0, 1.0); j++; // center for input 8
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for input 8
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; //radius for 
        //BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // radius for cosine term
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // weight for output 1
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // weight for output 2
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // weight for output 3
        BORG_Problem_set_bounds(problem, j, 0.0, 1.0); j++; // weight for output 4
    }
    
    BORG_Problem_set_bounds(problem, j, 0.0, 2*PI); j++; // phi-- phase shift
    //BORG_Problem_set_bounds(problem, j, 0.0, 2*PI); j++; // phi-- phase shift
 
    BORG_Problem_set_epsilon(problem, 0, 100); // to be fixed
    BORG_Problem_set_epsilon(problem, 1, 1000); // to be fixed
    BORG_Problem_set_epsilon(problem, 2, 0.1); // to be fixed
    BORG_Problem_set_epsilon(problem, 3, 1); // to be fixed
    //BORG_Problem_set_epsilon(problem, 5, 1); // to be fixed
    // This is set up to run only one seed at a time.
    char outputFilename[256];
    char runtime[256];
    FILE* outputFile = NULL;
    sprintf(outputFilename, "../sets/CRB_S%d.set", seed); // output path (make sure this exists)
    sprintf(runtime, "../runtime/CRB_S%d_M%%d.runtime", seed); // runtime path (make sure this exists)

    BORG_Algorithm_output_runtime(runtime);

    int rank; // different seed on each processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    BORG_Random_seed(37*seed*(rank+1));
    BORG_Archive result = BORG_Algorithm_ms_run(problem); // this actually runs the optimization

    // If this is the master node, print out the final archive
    if (result != NULL) {
        outputFile = fopen(outputFilename, "w");
        if (!outputFile) {
            BORG_Debug("Unable to open final output file\n");
        }
        BORG_Archive_print(result, outputFile);
        BORG_Archive_destroy(result);
        fclose(outputFile);
    }

    BORG_Algorithm_ms_shutdown();
    BORG_Problem_destroy(problem);


    return EXIT_SUCCESS;

}
  
