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

//Structures initialisation
Node *nodes;
Prob *pr;
Res *res;
int n, lb, ub, lf, uf, D;
double lc, uc;

int main(int arg, char * argv[])
{
    // Problem file
    char * instance_path = "./instances/g2_n1000_beta10_9";

    // Read instance file and initialize the problem
    nodes = read_file(instance_path, &n, &D, &lb, &ub, &lc, &uc, &lf, &uf);
    pr    = generate_prob(D, n, lb, ub, lc, uc, lf, uf);
    res   = create_res();

    // Solve the problem
    solveSSFCTP(nodes, pr, res);

    // Optimal value
    printf("opt_val = %lf\n", res->opt_val);

    // Clean memory
    free(nodes);
    free(pr);
    free(res);
}