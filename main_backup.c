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
int do_posCostsRed = 0; // Minimal unit cost c_j centered to 0             if true
int do_sort        = 1; // Phase 0 executed                                if true
int do_ZLP_ZG      = 1; // Phase 1 executed                                if true
int do_P3          = 1; // Phases 2 and 3 executed                         if true
int do_dom         = 1; // Phase 6 exectured                               if true
int do_LP          = 1; // Phase 7 executed                                if true
int do_filter      = 1; // Phase 9 executed                                if true
int do_exact       = 1; // Phases 10 executed                              if true
int computeSol     = 1; // Computes the optimal solution                   if true, just the optimal value otherwise

int verbose        = 0; // Intermediate results printed on standard output if true
int outputSol      = 0; // Prints the optimal solution on standard output  if true
int N              = 10; // Number of instances generated and solved
enum GenType gen = PisingerWeakDiff; // GenType (must be in {simple, Klose1, Klose2, correlated, PisingerStrong, PisingerWeak, File, PisingerWeakDiff, simpleDiff})

//Output options
int produceFileKlose   = 1; //if activated, a file containing the instance is saved and can be solved using Klose's algorithm later
int produceOutputJulia = 0; //if activated, priduces standard output to solve the problem in Julia

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

    //----------------------------------------------------- Problem parameters --------------------------------------------------
    //Random number generator seed
    int instance = 21;  // simpleDiff(M=0.5)->1, (M=1)->2, (M=5)->3, (M=10)->4, (M=25)->5 // PisingerWeakDiff(500,5)->6, (1000,5)->7, (5000,5)->8, (10000,5)->9, (500,10)->10, (1000,10)->11, (5000,10)->12, (10000,10)->13, (500,100)->14, (1000,100)->15, (5000,100)->16, (10000,100)->17, (500,1000)->18, (1000,1000)->19, (5000,1000)->20, (10000,1000)->21,
    int seed = initSeed((int)(123456789*instance)); // seed will be randomly generated if the seed given in argument is 0 

    //Parameters generator_simple
    if(gen == simple || gen == simpleDiff){
        double M = 25;         // {0.5, 1, 5, 10, 25}
        n =  (int)(1000*M);
        lb = (int)(500*M);
        ub = (int)(1000*M);
        lc = 8;
        uc = 12;
        lf = (int)(5000*M);
        uf = (int)(10000*M);
        //D = 500000*M*M; //Used in simple, but unused for simpleDiff
    }

    //Parameters for Klose1 instances
    else if(gen == Klose1){
        n = 10000;
        Br = 50;
        Fr = 1.00;
    }
        
    //Parameters for Klose2 instances
    else if(gen == Klose2){
        n = 50;
    }

    //Parameters for Pisinger instances
    else if(gen == correlated || gen == PisingerStrong || gen == PisingerWeak || gen == PisingerWeakDiff){
        n = 10000;        //{500, 1000, 5000, 10000}
        lb = 5000;      //{(5000,10000)}//
        ub = 10000;     /////////////////////////////////////////////
        beta = 1000;    //{5, 10, 100, 1000}
        theta = 0.3;    //{(0.3, 0.5)}//
        omega = 0.5;    ////////////////////////////////////////////////////
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
        nodes = generateNodes(nodes, &D, &lb, &ub, &lc, &uc, &lf, &uf, &n, &Br, &Fr, &seed, &i, &j, &N, &beta, &theta, &omega, gen);
        //printNodes(nodes,n);

	    //------------------ Problem generation ------------------
        pr = generate_prob(D, n, lb, ub, lc, uc, lf, uf);
        //printProblem(pr);

        //-------------------- Output options --------------------

        if(produceFileKlose) //Save the instance in a file that is usable by Klose's algorithms
	    	outputFileKlose(nodes, n, pr->D, i);   //File output to solve the problem with Klose's functions 

	    if(produceOutputJulia) //Standard output to solve the problem in Julia
	    	outputJulia(nodes, n);

	    //------------------ Problem resolution ------------------
        res = create_res();
        solveSSFCTP_detailed(nodes, pr, res, do_posCostsRed, do_sort, do_ZLP_ZG, do_P3, do_dom, do_LP, do_filter, do_exact, computeSol, verbose, outputSol);      
        updateGlobalRes(gRes, res);

        //--------------------- Clean memory ---------------------
        cleanMemory(res, P3L, pr, nodes);
    }
    //--------------------- Global stats --------------------
    printGlobalRes(gRes, computeSol);

    //--------------------- Clean memory --------------------
    free(gRes);
}



/*
 ===========================================================
|| ------------------- simpleDiff ------------------------ ||
||Instances non-corrélées (inspirées de Klose)             ||
||b_j is an uniform integer in [lb, ub]                    ||
||c_j is an uniform double in (lc, uc)                     ||
||f_j is an uniform integer in [lf, uf]                    ||
||                                                         ||
||n = 1000*M;                                              ||
||lb = 500*M;                                              ||
||ub = 1000*M;                                             ||
||lc = 8;                                                  ||
||uc = 12;                                                 ||
||lf = 5000*M;                                             ||
||uf = 10000*M;                                            ||
||Di = i/(N+1) sum(b_j) pour l'instance i \in {1,2,...,N}  ||
||                                                         ||
||paramétrisations testées:                                ||
||m              \in {0.5, 1, 5, 10}                       ||
 ===========================================================

 ===============================================================================
|| ----------------------------- PisingerWeakDiff ---------------------------- ||
||Instances corrélées (inspirées de Pisinger)                                  ||
||b_j is an uniform integer in [lb, ub]                                        ||
||u_j = b_j + alpha_j * ub/beta, where alpha_j is a uniform double on (-1,1),  ||
||beta is a constant (10/20/100/1000).                                         ||
||The bigger is beta, the more correlated are the instances                    ||
||                                                                             ||
||f_j = is a uniform integer in [theta*u_j, omega*u_j]                         ||
||c_j = (u_j-f_j)/b_j                                                          ||
||Di  = i/(N+1) sum(b_j) pour l'instance i \in {1,2,...,N}                     ||
||                                                                             ||
||paramétrisations testées:                                                    ||
||n              \in {500, 1000, 5000, 10000}                                  ||
||beta           \in {5, 10, 100, 1000}                                        ||
||(lb,ub)        \in {(5000,10000)}                                            ||
||(theta,omega)  \in {(0.3, 0.5)}                                              ||
 ===============================================================================

*/