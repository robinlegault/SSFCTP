#include "util.h"
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>

#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS 1
    #define restrict __restrict
#endif


#include <stdint.h>
#include <errno.h>

/*
POSIX getline replacement for non-POSIX systems (like Windows)
Differences:
    - the function returns int64_t instead of ssize_t
    - does not accept NUL characters in the input file
Warnings:
    - the function sets EINVAL, ENOMEM, EOVERFLOW in case of errors. The above are not defined by ISO C17,
    but are supported by other C compilers like MSVC
*/
int64_t my_getline(char **restrict line, size_t *restrict len, FILE *restrict fp) {
    // Check if either line, len or fp are NULL pointers
    if(line == NULL || len == NULL || fp == NULL) {
        errno = EINVAL;
        return -1;
    }
    
    // Use a chunk array of 128 bytes as parameter for fgets
    char chunk[128];

    // Allocate a block of memory for *line if it is NULL or smaller than the chunk array
    if(*line == NULL || *len < sizeof(chunk)) {
        *len = sizeof(chunk);
        if((*line = malloc(*len)) == NULL) {
            errno = ENOMEM;
            return -1;
        }
    }

    // "Empty" the string
    (*line)[0] = '\0';

    while(fgets(chunk, sizeof(chunk), fp) != NULL) {
        // Resize the line buffer if necessary
        size_t len_used = strlen(*line);
        size_t chunk_used = strlen(chunk);

        if(*len - len_used < chunk_used) {
            // Check for overflow
            if(*len > SIZE_MAX / 2) {
                errno = EOVERFLOW;
                return -1;
            } else {
                *len *= 2;
            }
            
            if((*line = realloc(*line, *len)) == NULL) {
                errno = ENOMEM;
                return -1;
            }
        }

        // Copy the chunk to the end of the line buffer
        memcpy(*line + len_used, chunk, chunk_used);
        len_used += chunk_used;
        (*line)[len_used] = '\0';

        // Check if *line contains '\n', if yes, return the *line length
        if((*line)[len_used - 1] == '\n') {
            return len_used;
        }
    }

    return -1;
}

//============================================= Comparators and sorting ======================================================

// Return 1 if e_1 < e_2, where e_j is the linearized cost of each unit on node j
// n1 and n2 are pointers to nodes
static int node_comparator(const void *n1, const void *n2)
{
  const Node *p1 = (Node *)n1;
  const Node *p2 = (Node *)n2;
  if (p1->e < p2->e)
    return -1;
  else if (p1->e > p2->e)
    return +1;
  else
    return 0;
}

void sort_nodes(int n, Node *nodes){
  qsort(nodes, n, sizeof(Node), node_comparator);            //Sorting the nodes in non-decreasing order of e_j
        for (int i = 0; i < n; ++i)                                //Node index fits the sorting order
            nodes[i].i = i;
}

//================================================= Nodes generation ==========================================================
// Generate a simple instance with n nodes. 
// For each node j, the capacity b_j, the unit cost c_j and the fixed cost f_j are defined as follows:
// b_j is an uniform integer in [lb, ub]
// c_j is an uniform double in (lc, uc)
// f_j is an uniform integer in [lf, uf]
//
// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_simple(int n, int lb, int ub, double lc, double uc, int lf, int uf, int seed)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(lb,ub);
    nodes[i].c = randDouble(lc,uc);
    nodes[i].f = randInt(lf,uf);
    nodes[i].u = nodes[i].f + nodes[i].b * nodes[i].c;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;
  }

  sort_nodes(n, nodes);
  return nodes;
}


// Same as generator_simple, but Demand for instance instanceIdx out of N is given by D= floor(instanceIdx/N+1 * sumBj),
// for instanceIdx in {1,2,3,...,N}
Node* generator_simple_different_demand(int* D, int n, int lb, int ub, double lc, double uc, int lf, int uf, int instanceIdx, int N, int seed)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));
  int sumBj = 0;

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(lb,ub);
    sumBj += nodes[i].b;
    nodes[i].c = randDouble(lc,uc);
    nodes[i].f = randInt(lf,uf);
    nodes[i].u = nodes[i].f + nodes[i].b * nodes[i].c;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;
  }

  *D = floor((double)(instanceIdx)/(N+1)*sumBj);
  sort_nodes(n, nodes);
  return nodes;
}


// Generate an instance with n nodes using Klose's group 1 problems generation method, using the two following parameters:
//     - The B-ratio, Br = 100*D/Sum{b_j}, with D = 100000. Br must be between in (0,100), with Br = 50 to have approximately 50% of the nodes in the optimal solution
//     - The F-ratio, Fr = mean(f_j)/(mean(c_j)*mean(b_j)). Fr must be between in [0,+inf), with Fr = 1 to have approximately half of the total cost from fixed cost and half from partial costs
//
// For each node j, the capacity b_j, the unit cost c_j and the fixed cost f_j are defined as follows:
// b_j is an uniform integer in [7500, 12500]
// c_j is an uniform double in (8, 12)
// f_j is an uniform integer in [75000, 125000]
// The capacities and the fixed costs are then scaled to respect the given B-ratio and F-ratio
//
// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_klose1(int* D, int* lb, int* ub, double* lc, double* uc, int* lf, 
  int* uf, int n, double Br, double Fr, int seed)
{
  *D  = 100000; 
  *lb = 7500; 
  *ub = 12500;
  *lc = 8;
  *uc = 12;
  *lf = 75000;
  *uf = 125000;
  int sumB = 0;
  int sumF = 0;
  double   sumC = 0;

  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  int i;
  for(i=0; i<n; i++){
    nodes[i].b = randInt(*lb,*ub);
    sumB += nodes[i].b;
    nodes[i].c = randDouble(*lc,*uc);
    sumC += nodes[i].c;
    nodes[i].f = randInt(*lf,*uf);
    sumF += nodes[i].f;
  }

  double targetSumB = 100*(*D)/Br;
  double mult = targetSumB/sumB;
  sumB = 0;

  for (int i = 0; i < n; ++i)
  {
    nodes[i].b = round(nodes[i].b * mult);
    sumB += nodes[i].b ;
  }

  double meanB = sumB/n;
  double meanC = sumC/n;
  double meanF = sumF/n;
  double targetMeanF = meanB * meanC * Fr;

  mult = targetMeanF/meanF;

  for (int i = 0; i < n; ++i)
  {
    nodes[i].f = round(nodes[i].f * mult);

    nodes[i].u = nodes[i].f + nodes[i].b * nodes[i].c;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].i = i;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;
  }

  sort_nodes(n, nodes);
  return nodes;
}

// Generate an instance with n nodes using Klose's group 2 problems generation method.
//For each node j, the capacity b_j, the total cost u_j, the unit cost c_j and the fixed cost f_j are defined as follows:

//b_j is an uniform integer in [10, 100]. 
//u_j is a uniform integer in [b_j, 2*b_j]
//f_j = phi_j*u_j, where phi_j is a uniform double in (0.75, 1)
//c_j = (u_j-f_j)/b_j is set so the total cost is really u_j on node j
//The demand D is fixed to 0.5*Sum{b_j} (B-ratio is 50)

//The linearized cost e_j = u_j/b_j on each node lies between 1 and 2, with an expected value of 1.5

// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_klose2(int* D, int* lb, int* ub, double* lc, double* uc, int* lf, 
  int* uf, int n, int seed)
{
  *lb = 10; 
  *ub = 100;
  int sumB = 0;
  double totCost;
  double phi_i;
  int beta_i;

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;

  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  int i;
  for(i=0; i<n; i++){
    nodes[i].b = randInt(*lb,*ub);
    sumB += nodes[i].b;
    beta_i = randInt(0,nodes[i].b);
    totCost = nodes[i].b + beta_i;
    phi_i = randDouble(0.75,1);
    nodes[i].f = round(phi_i*totCost);
    nodes[i].c = (totCost - nodes[i].f)/nodes[i].b;
    nodes[i].u = nodes[i].f + nodes[i].b*nodes[i].c;
    nodes[i].e = nodes[i].u/nodes[i].b;

    nodes[i].i = i;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;

    if(nodes[i].c < minC){
      minC = nodes[i].c;
    }
    if(nodes[i].c > maxC){
      maxC = nodes[i].c;
    }
    if(nodes[i].f < minF){
      minF = nodes[i].f;
    }
    if(nodes[i].f > maxF){
      maxF = nodes[i].f;
    }
  }

  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;
  *D = round(0.5*sumB);
  sort_nodes(n, nodes);
  return nodes;
}


// Correlated instance (inspired by the correlated instance from "Where are the hard knapsack problems? - Pisinger 2005")
// For each node j, the capacity b_j, the unit cost c_j, the fixed cost f_j and the total cost u_j are defined as follows:
// b_j is an uniform integer in [lb, ub]
// u_j = b_j*10 + ub
// f_j = is a uniform integer in [0.3*u_j, 0.5*u_j]
// c_j = (u_j-f_j)/b_j
//
// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_correlated(int* D, int* lb, int* ub, double* lc, double* uc, int* lf, 
  int* uf, int n, int seed)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  *lb = 750; 
  *ub = 1250;
  *D = 500*n;

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(*lb,*ub);
    nodes[i].u = nodes[i].b*10 + randDouble(0.9,1)*(*ub);
    nodes[i].f = randInt(ceil(0.3*nodes[i].u),ceil(0.5*nodes[i].u));
    nodes[i].c = (nodes[i].u - nodes[i].f)/nodes[i].b;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;

    if(nodes[i].c < minC){
      minC = nodes[i].c;
    }
    if(nodes[i].c > maxC){
      maxC = nodes[i].c;
    }
    if(nodes[i].f < minF){
      minF = nodes[i].f;
    }
    if(nodes[i].f > maxF){
      maxF = nodes[i].f;
    }
  }

  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  sort_nodes(n, nodes);
  return nodes;
}

// Strongly correlated instance (inspired by the strongly correlated instance from "Where are the hard knapsack problems? - Pisinger 2005")
// For each node j, the capacity b_j, the unit cost c_j, the fixed cost f_j and the total cost u_j are defined as follows:
// b_j is an uniform integer in [1, ub]
// u_j = b_j + ub/10
// f_j = is a uniform integer in [0.3*u_j, 0.5*u_j]
// c_j = (u_j-f_j)/b_j
// D = (sum{b_j})/3
//
// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_Pisinger_strongly_correlated(int* D, int lb, int ub, double* lc, double* uc, int* lf, 
  int* uf, int n, int seed)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  *D = round((double)(ub+lb)/3)*n;

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(lb,ub);
    nodes[i].u = nodes[i].b + 0.1*ub;
    nodes[i].f = randInt(ceil(0.3*nodes[i].u),ceil(0.5*nodes[i].u));
    nodes[i].c = (nodes[i].u - nodes[i].f)/nodes[i].b;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;

    if(nodes[i].c < minC){
      minC = nodes[i].c;
    }
    if(nodes[i].c > maxC){
      maxC = nodes[i].c;
    }
    if(nodes[i].f < minF){
      minF = nodes[i].f;
    }
    if(nodes[i].f > maxF){
      maxF = nodes[i].f;
    }
  }

  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  sort_nodes(n, nodes);
  return nodes;
}

// Weakly correlated instance (inspired by the weakly correlated instance from "Where are the hard knapsack problems? - Pisinger 2005")
// For each node j, the capacity b_j, the unit cost c_j, the fixed cost f_j and the total cost u_j are defined as follows:
// b_j is an uniform integer in [lb, ub]
// u_j = b_j + alpha_j * ub/beta, where alpha_j is a uniform double on (-1,1), beta is a constant (10/20/100/1000). The bigger is beta, the more correlated are the instances
// f_j = is a uniform integer in [0.3*u_j, 0.5*u_j]
// c_j = (u_j-f_j)/b_j
// D = (sum{b_j})/3
//
// Nodes are returned sorted by non-decreasing order of e_j
Node* generator_Pisinger_weakly_correlated(int* D, int lb, int ub, double* lc, double* uc, int* lf, 
  int* uf, int n, int seed)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  *D = round((double)(ub+lb)/3)*n;

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(lb,ub);
    nodes[i].u = nodes[i].b + randDouble(-1,1)*0.1*ub;
    nodes[i].f = randInt(ceil(0.3*nodes[i].u),ceil(0.5*nodes[i].u));
    nodes[i].c = (nodes[i].u - nodes[i].f)/nodes[i].b;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;

    if(nodes[i].c < minC){
      minC = nodes[i].c;
    }
    if(nodes[i].c > maxC){
      maxC = nodes[i].c;
    }
    if(nodes[i].f < minF){
      minF = nodes[i].f;
    }
    if(nodes[i].f > maxF){
      maxF = nodes[i].f;
    }
  }

  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  sort_nodes(n, nodes);
  return nodes;
}


// Same as generator_Pisinger_weakly_correlated, but Demand for instance instanceIdx out of N is given by D= floor(instanceIdx/N+1 * sumBj),
// for instanceIdx in {1,2,3,...,N}. Also, parameters beta, theta and omega can be specified
Node* generator_Pisinger_weakly_correlated_different_demand(int* D, int lb, int ub, double* lc, double* uc, int* lf, 
  int* uf, int n, int seed, int instanceIdx, int N, double beta, double theta, double omega)
{
  srand( (unsigned) seed);
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;
  int sumBj = 0;

  int i;
  for(i=0; i<n; i++){
    nodes[i].i = i;
    nodes[i].b = randInt(lb,ub);
    sumBj += nodes[i].b;
    nodes[i].u = nodes[i].b + randDouble(-1,1)/beta*ub;
    nodes[i].f = randInt(ceil(theta*nodes[i].u),ceil(omega*nodes[i].u));
    nodes[i].c = (nodes[i].u - nodes[i].f)/nodes[i].b;
    nodes[i].e = nodes[i].u / nodes[i].b;
    nodes[i].LBx = 1;
    nodes[i].UBx = nodes[i].b-1;
    nodes[i].LB = 0;
    nodes[i].poss_p = 1;
    nodes[i].poss_1 = 1;
    nodes[i].poss_0 = 1;
    nodes[i].x = 0;
    nodes[i].y = 0;
    nodes[i].z = 0;

    if(nodes[i].c < minC){
      minC = nodes[i].c;
    }
    if(nodes[i].c > maxC){
      maxC = nodes[i].c;
    }
    if(nodes[i].f < minF){
      minF = nodes[i].f;
    }
    if(nodes[i].f > maxF){
      maxF = nodes[i].f;
    }
  }

  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  sort_nodes(n, nodes);

  *D = floor((double)(instanceIdx)/(N+1)*sumBj);
  return nodes;
}

// Generate an instance using file 1_j_i.txt in repertory "instances"
Node* generator_file(int j, int i, int *D, int* lb, int* ub, double* lc, double* uc, int* lf, int* uf, int* n)
{
  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  //instances//
  char pname[1024]= "./instances/1_"; 
  char str[1024];
  sprintf(str, "%d", j);
  strcat( pname, str );
  strcat( pname, "_" );
  sprintf(str, "%d", i);
  strcat( pname, str );
  strcat( pname, ".txt");
  printf("nom fichier = %s\n", pname);

  //snft10x100//
  //char pname[1024]= "./snft10x100/test10x1000_"; 
  //char str[1024];
  //sprintf(str, "%d", j);
  //strcat( pname, str );
  //strcat( pname, "_" );
  //sprintf(str, "%d", i);
  //strcat( pname, str );
  //strcat( pname, ".txt");
  //printf("nom fichier = %s\n", pname);

  fp = fopen(pname, "r");
  if (fp == NULL)
      exit(EXIT_FAILURE);

  read = my_getline(&line, &len, fp);//TODO
  char * token = strtok(line, " ");
  sscanf(token, "%d", n); 
  token = strtok(NULL, " ");
  sscanf(token, "%d", D); 

  double minC = DBL_MAX;
  double maxC  = 0.0;
  int minF  = INT_MAX;
  int maxF  = 0.0;
  int minB  = INT_MAX;
  int maxB  = 0.0;

  //Node *nodes  = (Node*)malloc((*n) * sizeof(Node)); 
  Node *nodes = (Node*)calloc((*n), sizeof(Node));

  for (int idx = 0; idx < *n; ++idx)
  {
    (read = my_getline(&line, &len, fp));//TODO
      double c;
      int f;
      int b;
      sscanf(line, "%lf %d %d", &c, &f, &b);
      double u = (c*b+f);
      double e = c+f/b;

      nodes[idx].i = idx;
      nodes[idx].b = b;
      nodes[idx].u = u;
      nodes[idx].f = f;
      nodes[idx].c = c;
      nodes[idx].e = e;
      
      //printf("-------------------------------------\n");
      //printf("e[%d]     = %f\n",idx,e );
      //printf("c[%d]     = %f\n",idx,c );
      //printf("nodes[%d].e = %f\n",idx,nodes[idx].e );
      //printf("nodes[%d].c = %f\n",idx,nodes[idx].c );
      //printf("-------------------------------------\n");
      
      nodes[idx].LBx = 1;
      nodes[idx].UBx = b-1;
      nodes[idx].LB = 0;
      nodes[idx].poss_p = 1;
      nodes[idx].poss_1 = 1;
      nodes[idx].poss_0 = 1;
      nodes[idx].x = 0;
      nodes[idx].y = 0;
      nodes[idx].z = 0;

      if(nodes[idx].c < minC){
        minC = nodes[idx].c;
      }
      if(nodes[idx].c > maxC){
        maxC = nodes[idx].c;
      }
      if(nodes[idx].f < minF){
        minF = nodes[idx].f;
      }
      if(nodes[idx].f > maxF){
        maxF = nodes[idx].f;
      }
      if(nodes[idx].b < minB){
        minB = nodes[idx].b;
      }
      if(nodes[idx].b > maxB){
        maxB = nodes[idx].b;
      }
  }

  fclose(fp);
  if (line)
      free(line);
  
  *lb = minB;
  *ub = maxB;
  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  sort_nodes(*n, nodes);
  return nodes;
  //return NULL;
}


// Generate 
void generator_FCTP(int m, int n, int la, int ua, int lb, int ub, int lc, int uc, int lf, int uf, int seed)
{
  FILE * outfile;
  char pname[1024]= "./instance_FCTP/"; 
  char str[1024];
  strcat( pname, "ran" );
  sprintf(str, "%d", m);
  strcat( pname, str ); 
  strcat( pname, "x" );
  sprintf(str, "%d", n);
  strcat( pname, str ); 
  strcat( pname, ".dat" );
  //printf("Generating file %s\n",pname);
  outfile= fopen( pname, "wt" );
  fprintf(outfile, "%d\n", m );
  fprintf(outfile, "%d\n\n\n", n );

  int sum_ai =0;
  int sum_bj =0;
  int ai,bj,cij,fij;

  int a[m];
  int b[n];

  srand( (unsigned) seed);

  for (int j = 0; j < n; ++j)
  {
    b[j] = randInt(lb,ub);
    sum_bj += b[j];
  }

  for (int i = 0; i < m-1; ++i)
  {
    a[i] = randInt(la,ua);
    sum_ai += a[i];
  }

  a[m-1] = sum_bj - sum_ai;

  for (int i = 0; i < m; ++i)
  {
    fprintf(outfile, "%d\n", a[i] );
  }
  fprintf(outfile, "\n\n" );

  for (int j = 0; j < n; ++j)
  {
    fprintf(outfile, "%d\n", b[j] );
  }
  fprintf(outfile, "\n\n" );

  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j){
      fprintf(outfile, "%d %d %d %d\n", i, j, randInt(lc,uc), randInt(lf,uf));
    }
  }

  fclose( outfile );
}

//================================================= Problem generation ==========================================================
// Generate a problem (prob struct) corresponding to the given parameters
Prob* generate_prob(int D, int n, int lb, int ub, double lc, double uc, int lf, int uf){
  Prob *pr = (Prob*)calloc(1, sizeof(Prob));
  pr->D = D;
  pr->n = n;
  pr->lb = lb;
  pr->ub = ub;
  pr->lc = lc;
  pr->uc = uc;
  pr->lf = lf;
  pr->uf = uf;
  pr->UB = DBL_MAX;
  pr->LB = 0;
  pr->Cmin = lc;
  pr->Cmax = uc;
  pr->s  = 0;
  pr->xs = 0;
  pr->LB_linear = 0;
  pr->nb_poss_p = n;
  pr->nb_poss_1 = n;
  pr->nb_poss_0 = n;
  pr->cmin_closed_interval=1;
  pr->cmax_closed_interval=1;
  pr->incumbent = (int*)calloc(n, sizeof(int));
  return pr;
}

Node* generate_nodes(Node* nodes, int *D, int *lb, int *ub, double *lc, double *uc, int *lf, int *uf, int *n, double *Br, double *Fr, int *seed, int *i, int *j, int *N, double *beta, double *theta, double *omega, enum GenType gen){
  if(gen == simple){
      nodes = generator_simple(*n, *lb, *ub, *lc, *uc, *lf, *uf, *seed);    
  }
  else if(gen == Klose1){                            
      nodes = generator_klose1(D, lb, ub, lc, uc, lf, uf, *n, *Br, *Fr, *seed);
  }
  else if(gen == Klose2){            
      nodes = generator_klose2(D, lb, ub, lc, uc, lf, uf, *n, *seed);
  }
  else if(gen == correlated){                    
      nodes = generator_correlated(D, lb, ub, lc, uc, lf, uf, *n, *seed);
  }
  else if(gen == PisingerStrong){             
      nodes = generator_Pisinger_strongly_correlated(D, *lb, *ub, lc, uc, lf, uf, *n, *seed);
  }
  else if(gen == PisingerWeak){   
      nodes = generator_Pisinger_weakly_correlated(D, *lb, *ub, lc, uc, lf, uf, *n, *seed);
  }
  else if(gen == File){     
      nodes = generator_file(*j, *i, D, lb, ub, lc, uc, lf, uf, n); // file 1_j_i.txt
  }
  else if(gen == PisingerWeakDiff){
      nodes = generator_Pisinger_weakly_correlated_different_demand(D, *lb, *ub, lc, uc, lf, uf, *n, *seed, *i+1, *N, *beta, *theta, *omega);
  }
  else if(gen == simpleDiff){
      nodes = generator_simple_different_demand(D, *n, *lb, *ub, *lc, *uc, *lf, *uf, *i+1, *N, *seed);
  }
  else{
      printf("Error: invalid generator type \n");
      return NULL;
  }
  return nodes;
}


