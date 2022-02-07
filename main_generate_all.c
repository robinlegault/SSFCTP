#include "combo.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "generator.h"
#include "combo.h"
#include "ssfctp.h"
#include <time.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <float.h>
#include <limits.h>

//User-defined options
int do_posCostsRed          = 1; // Minimal unit cost c_j centered to 0                                    if true
int do_sort                 = 1; // Sorting executed                                                       if true
int do_ZLP_ZG               = 1; // H1 executed                                                            if true
int do_P3                   = 1; // H2 executed                                                            if true
int do_dom                  = 1; // H3 exectured                                                           if true
int do_LP                   = 1; // H4 executed                                                            if true
int do_filter               = 1; // F executed                                                             if true
int do_exact                = 1; // E executed                                                             if true
int computeSol              = 1; // Computes the optimal solution                                          if true, just the optimal value otherwise
int alternative_lambda_init = 0; // Use max_{j in {1,2,...,s}}{c_j} as initial value of lambda in P3Search if true, and min{e_s, max_{j in {1,2,...,n}}{c_j}} otherwise

int verbose        = 1; // Intermediate results printed on standard output if verbose=2, only optimal value and cpu time if verbose=1, nothing if verbose=0
int outputSol      = 0; // Prints the optimal solution on standard output  if true
int N              = 10; // Number of instances generated and solved
//enum GenType gen = simpleDiff; // GenType (must be in {simple, Klose1, Klose2, correlated, PisingerStrong, PisingerWeak, File, PisingerWeakDiff, simpleDiff})
char gen_str[50] = "...";

//Article notation
// Group 1 (Uncorrelated data instances) -> simpleDiff
// Group 2 (Correlated data instances)   -> PisingerWeak
// Group 3 (Uncorrelated data instances) -> Klose1
// Group 4 (Uncorrelated data instances) -> Klose2

//Output options
int produceFileKlose   = 0; //if activated, a file containing the instance is saved and can be solved using Klose's algorithm later
int produceFileGurobi  = 0; //if activated, a file containing the instance is saved and can be solved using Gurobi later
int produceOutputJulia = 0; //if activated, priduces standard output to solve the problem in Julia 
int produceFileFinal = 1; //if activated, a file containing the instance is saved

//Structures initialisation
P3struct *P3L;
Node *nodes;
Prob *pr;
Res *res;
GlobalRes *gRes;

//Variables
int n=0, lb=0, ub=0, lf=0, uf=0, D=0, j=0;
double lc=0, uc=0, fixedObjectiveCost=0, Br=0, Fr=0, beta=0, theta=0, omega=0;

int main(int arg, char * argv[])
{
    enum GenType list_gen[77] = {simpleDiff, simpleDiff, simpleDiff, simpleDiff, simpleDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, PisingerWeakDiff, Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose1,Klose2,Klose2,Klose2,Klose2};
    int list_instance[77]= {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,251,252,253,254,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100};
    double list_M[77] = {0.5, 1, 5, 10, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int list_n[77] = {0,0,0,0,0, 500, 1000, 5000, 10000, 500, 1000, 5000, 10000, 500, 1000, 5000, 10000, 500, 1000, 5000, 10000, 25000, 25000, 25000, 25000,500,500,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,10000,500,1000,5000,10000};
    double list_beta[77] = {0,0,0,0,0,5,5,5,5,10,10,10,10,100,100,100,100,1000,1000,1000,1000,5,10,100,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double list_Br[77] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,10,10,10,25,25,25,50,50,50,5,5,5,10,10,10,25,25,25,50,50,50,5,5,5,10,10,10,25,25,25,50,50,50,5,5,5,10,10,10,25,25,25,50,50,50,0,0,0,0};
    double list_Fr[77] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0.3,0.6,1.0,0,0,0,0};

    for (int it = 0; it < 77; ++it) // Generate and solve N instances
    {
        enum GenType gen = list_gen[it];
        int instance = list_instance[it];
        double M = list_M[it];
        n = list_n[it];
        beta = list_beta[it];
        Br = list_Br[it];
        Fr = list_Fr[it];
        char instance_name[1024] = "";

        //----------------------------------------------------- Problem parameters --------------------------------------------------
        //Random number generator seed
        //int instance = 1;  // simpleDiff(M=0.5)->1, (M=1)->2, (M=5)->3, (M=10)->4, (M=25)->5 // PisingerWeakDiff(500,5)->6, (1000,5)->7, (5000,5)->8, (10000,5)->9, (500,10)->10, (1000,10)->11, (5000,10)->12, (10000,10)->13, (500,100)->14, (1000,100)->15, (5000,100)->16, (10000,100)->17, (500,1000)->18, (1000,1000)->19, (5000,1000)->20, (10000,1000)->21, (25000,5)->251, (25000,10)->252, (25000,100)->253, (25000,1000)->254
        int seed = initSeed((int)(123456789*instance)); // seed will be randomly generated if the seed given in argument is 0 // seed for section 7.3: instance = 1

        //Parameters generator_simple
        if(gen == simple || gen == simpleDiff){
            //double M =0.5;         // {0.5, 1, 5, 10, 25}
            printf("....................................................\n");
            printf("instance = %d, M = %lf\n", instance, M);
            printf("....................................................\n");
            n =  (int)(1000*M);
            lb = (int)(500*M);
            ub = (int)(1000*M);
            lc = 8;
            uc = 12;
            lf = (int)(5000*M);
            uf = (int)(10000*M);
            //D = 500000*M*M; //Used in simple, but unused for simpleDiff

            //instance name
            strcat(instance_name, "g1");
            char str[1024];
            sprintf(str, "_n%d", n);
            strcat(instance_name, str);
        }

        //Parameters for Klose1 instances
        else if(gen == Klose1){
            //n = 10000;    // 500 1000 5000 10000
            //Br = 50;    // 5 10 25 50
            //Fr = 1.0;  // 0.3 0.6 1.0
            printf("....................................................\n");
            printf("instance = %d, n = %d, Br = %lf, Fr = %lf\n", instance, n, Br, Fr);
            printf("....................................................\n");

            //instance name
            char str[1024];
            strcat(instance_name, "g3");
            sprintf(str, "_n%d", n);
            strcat(instance_name, str);
            sprintf(str, "_Br%d", (int)(Br));
            strcat(instance_name, str);
            sprintf(str, "_Fr%d", (int)(100*Fr));
            strcat(instance_name, str);
        }
            
        //Parameters for Klose2 instances
        else if(gen == Klose2){
            //n = 10000; // 500 1000 5000 10000
            printf("....................................................\n");
            printf("instance = %d, n = %d\n", instance, n);
            printf("....................................................\n");

            //instance name
            char str[1024];
            strcat(instance_name, "g4");
            sprintf(str, "_n%d", n);
            strcat(instance_name, str);
        }

        //Parameters for Pisinger instances
        else if(gen == correlated || gen == PisingerStrong || gen == PisingerWeak || gen == PisingerWeakDiff){
            //n     = 25000;      //{500, 1000, 5000, 10000, 25000}
            //beta  = 1000;        //{5, 10, 100, 1000}
            printf("....................................................\n");
            printf("instance = %d, n = %d, beta = %lf\n", instance, n, beta);
            printf("....................................................\n");
            lb    = 5000;     //{(5000,10000)}//
            ub    = 10000;    /////////////////////////////////////////////
            theta = 0.3;      //{(0.3, 0.5)}//
            omega = 0.5;      /////////////////////////////////////////////

            //instance name
            char str[1024];
            strcat(instance_name, "g2");
            sprintf(str, "_n%d", n);
            strcat(instance_name, str);
            sprintf(str, "_beta%d", (int)(beta));
            strcat(instance_name, str);
        }
        
        //Parameters generator_file
        else if(gen == File){
            int j = 504;
        }

        // Global results
        gRes = createGlobalRes();

        for (int i = 0; i < N; ++i) // Generate and solve N instances
        {

            //------------------- Nodes generation -------------------
            advanceSeed(&seed);
            nodes = generate_nodes(nodes, &D, &lb, &ub, &lc, &uc, &lf, &uf, &n, &Br, &Fr, &seed, &i, &j, &N, &beta, &theta, &omega, gen);
            //printNodes(nodes,n);

            //------------------ Problem generation ------------------
            pr = generate_prob(D, n, lb, ub, lc, uc, lf, uf);
            //printProblem(pr);

            //-------------------- Output options --------------------

            if(produceFileKlose) //Save the instance in a file that is usable by Klose's algorithms
                outputFileKlose(nodes, n, pr->D, i);   //File output to solve the problem with Klose's functions 

            if(produceOutputJulia) //Standard output to solve the problem in Julia
                outputJulia(nodes, n);

            if(produceFileGurobi) //Save the instance in a .mps file that is usable by Gurobi
                outputFileGurobi(nodes, n, pr->D, gen_str, i);   //File output to solve the problem with Klose's functions 

            if(produceFileFinal){ //Save the instance in a file (final format)
                outputFileFinal(nodes, n, pr->D, instance_name, i); //File output to solve the problem with Klose's functions 
            }

            //------------------ Problem resolution ------------------
            res = create_res();
            solveSSFCTP_detailed(nodes, pr, res, do_posCostsRed, do_sort, do_ZLP_ZG, do_P3, do_dom, do_LP, do_filter, do_exact, computeSol, verbose, outputSol, alternative_lambda_init);
            //solveSSFCTP(nodes, pr, res);  
            updateGlobalRes(gRes, res);

            //--------------------- Clean memory ---------------------
            cleanMemory(res, P3L, pr, nodes);
        }
        //--------------------- Global stats --------------------
        printGlobalRes(gRes, computeSol);

        //--------------------- Clean memory --------------------
        free(gRes);
    }
}