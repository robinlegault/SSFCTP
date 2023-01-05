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
enum GenType gen = PisingerWeakDiff; // GenType (must be in {simple, Klose1, Klose2, correlated, PisingerStrong, PisingerWeak, File, PisingerWeakDiff, simpleDiff})
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

//Structures initialisation
P3struct *P3L;
Node *nodes;
Prob *pr;
Res *res;
GlobalRes *gRes;

//Variables
int n=0, lb=0, ub=0, lf=0, uf=0, D=0, j=0;
double lc=0, uc=0, fixedObjectiveCost=0, Br=0, Fr=0, beta=0, theta=0, omega=0;

int n2[5]={500, 1000, 5000, 10000, 25000};
int beta2[4]={5, 10, 100, 1000};

int main(int arg, char * argv[])
{
    
    for (int idx_n = 0; idx_n < 5; idx_n++) // 
    {
        for (int idx_beta = 0; idx_beta < 4; idx_beta++) // 
        {
            // Results for the 10 instances of group 2 with n=n2[idx_n] and beta=beta2[idx_beta]
            gRes = createGlobalRes();
            printf("\n\nResults n=%d, beta=%d\n",n2[idx_n],beta2[idx_beta]);
            for (int j = 0; j < 10; ++j) // 
            {
                // Problem file
                char instance_path[1024]= "./instances/g2_n"; 
                char str[1024];
                sprintf(str, "%d", n2[idx_n]);
                strcat( instance_path, str ); 
                strcat( instance_path, "_beta" ); 
                sprintf(str, "%d", beta2[idx_beta]);
                strcat( instance_path, str ); 
                strcat( instance_path, "_" ); 
                sprintf(str, "%d", j);
                strcat( instance_path, str ); 
                printf("%s: ",instance_path);

                // Read instance file and initialize the problem
                nodes = read_file(instance_path, &n, &D, &lb, &ub, &lc, &uc, &lf, &uf);
                pr    = generate_prob(D, n, lb, ub, lc, uc, lf, uf);
                res   = create_res();

                //------------------ Problem resolution ------------------
                solveSSFCTP_detailed(nodes, pr, res, do_posCostsRed, do_sort, do_ZLP_ZG, do_P3, do_dom, do_LP, do_filter, do_exact, computeSol, verbose, outputSol, alternative_lambda_init);
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
}