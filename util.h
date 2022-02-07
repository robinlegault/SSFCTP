#ifndef UTIL_H
#define UTIL_H
#include <inttypes.h>

//================================================== Macros =========================================================
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
//================================================ Constants ========================================================

#define EPSILON 0.000000001   //used to eliminate rounding errors in comparisons
//#define PATH_MAX 4096

//============================================= Nb manipulations ====================================================
int ceil_log2(int x);

//======================================== Random number generation =================================================

enum GenType {simple, Klose1, Klose2, correlated, PisingerStrong, PisingerWeak, File, PisingerWeakDiff, simpleDiff};
double randU01();
int randInt(int min, int max);
double randDouble(double min, double max);

//============================================= Data structures =====================================================

//Contains relevant information on a single node from an instance of the ssfctp
typedef struct Node{
  int    i_init;     /* index (original) */
  int    i;          /* index (sorted) */
  int    f;      /* fixed cost */
  int    b;      /* offer */
  double c;      /* unit cost */
  double e;      /* linearized unit cost */
  double u;      /* total cost */
  int    LBx;    /* LB on (x | i is partial) */
  int    UBx;    /* UB on (x | i is partial) */
  double LB;     /* LB on (P3| i is partial) */
  int    poss_p; /* TRUE if the node may be partial (LBxp<UBxp), FALSE otherwise */
  int    poss_1; /* TRUE if the node may be completely used, FALSE otherwise */
  int    poss_0; /* TRUE if the node may be unused, FALSE otherwise */

  int  x; /* state of the node in the incumbent solution (may change later) */
  int       y;
  int       z;
} Node;

//Contains relevant information on an instance of the ssfctp
typedef struct Prob{
  int    D;  //Demand
  int    n;  //Number of nodes
  int    lb; //LB on b_j
  int    ub; //UB on b_j
  double lc; //LB on c_j
  double uc; //UB on c_j
  int    lf; //LB on f_j
  int    uf; //UB on f_j

  int    s;                    //Index of split node
  int    xs;                   //Usage of split node s in the linear relaxation's optimal solution 
  double LB_linear;            //Linear relaxation's optimal solution's objective cost 
  int    nb_poss_p;            //Number of possibly partial nodes
  int    nb_poss_1;            //Number of possibly completely used nodes
  int    nb_poss_0;            //Number of possibly unused nodes
  double LB;                   //LB on the optimal solution's objective function
  double UB;                   //UB on the optimal solution's objective function
  double Cmin;                 //LB on the partial node's unit cost c_j (if the optimal solution contains a partial node)
  double Cmax;                 //UB on the partial node's unit cost c_j (if the optimal solution contains a partial node)
  int    cmin_closed_interval; //1 if a node p with c_p=cmin can improve strictly the incumbent solution, 0 otherwise
  int    cmax_closed_interval; //1 if a node p with c_p=cmax can improve strictly the incumbent solution, 0 otherwise
  int*   incumbent;            //incumbent solution (x1*,x2*,...,xn*)
} Prob;

//Contains relevant information on the solution to a linear relaxation of the ssfctp
typedef struct LRstruct{
  double   cost;       /*Linear relaxation objective value*/
 //double   cp;       /*Maximal unit cost among nodes used in the relaxation solution*/
 // int      p;        /*Index of node p of unit cost cp*/
  int s;        /*Split item of the relaxation*/
  int xs;       /*Units used on s*/
} LRstruct;

typedef struct P3struct{
  double   lambda;   /* lambda */
  double   cp;       /* unit cost cp on p */
  int      p;        /* index of partial node p */
  double   P3lambda; /* objective value P3(lambda) */
  double   Zssfctp;  /* objective value P3(lambda) */
  int      surplus;
} P3struct;

typedef struct NodeLinkedList{
  struct Node* ptrNode;
  struct NodeLinkedList *prev;
  struct NodeLinkedList *next;
} NodeLinkedList;

typedef struct DominanceNodeLinkedList{
  struct Node* ptrNode; //pointer to the current dominant node on the current interval 
  double xMax;          //upper bound of the current interval (where this node is dominant)
  struct DominanceNodeLinkedList *next;//pointer to the dominant node on the next interval
} DominanceNodeLinkedList;

//=========================== Heap functions (to be used on an array of pointers) ===================================
int heap_insert(Node** nodesLR, Node* ptrNode, int* nbElementsLR);
Node* heap_pop(Node** nodesLR, int* nbElementsLR);
void printHeap(Node** nodesLR, int* nbElementsLR);
int partial_node_comparator(Node *p1, Node *p2);

//======================================== Build/clean functions for LBp2 ===========================================
void build_heap_LR(int n, int D, double *ZLR, int *s, Node *nodes, Node **nodesLR, int *nbElementsLR);
void build_linkedList_nodes(NodeLinkedList **firstNodeLL, NodeLinkedList **splitNodeLL, NodeLinkedList **lastNodeLL, Node *nodes, int s, int n);
void build_arraylist_partials(int n, Node * nodes, Node **nodesPart, int *nb_poss_part);
void util_free_LL(NodeLinkedList* firstNodeLL, NodeLinkedList* lastNodeLL);

//================================================ Output functions =================================================
void printNodesPart(Node **nodesPart, int nb_poss_part);

//=========================================== Detailed results container ============================================
//Contains detailed statistics on the execution of the algorithm
typedef struct Res{
  uint64_t cpu_sort_0;      /*cpu time spent                   in phase 0  (in microseconds)   */
  uint64_t cpu_ZLP_ZG_1;    /*cpu time spent                   in phase 1  (in microseconds)   */
  uint64_t cpu_P3_2_3_4_5;  /*cpu time spent             in phases 2 to 5  (in microseconds)   */
  uint64_t cpu_dom_6;       /*cpu time spent                   in phase 6  (in microseconds)   */
  uint64_t cpu_LP_7;        /*cpu time spent                   in phase 7  (in microseconds)   */
  uint64_t cpu_filtering_9; /*cpu time spent                   in phase 9  (in microseconds)   */
  uint64_t cpu_exact_10;    /*cpu time spent                   in phase 10 (in microseconds)   */
  uint64_t cpu_tot;         /*total cpu time spent, i.e. in phases 0 to 10 (in microseconds)   */
  int prOpt_ZLP_ZG_1;       /*equals 1 if the incumbent sol. is proved to be opt. after phase 1                */
  int isOpt_ZLP_ZG_1;       /*equals 1 if the incumbent sol. is optimal (not necessarily proven) after phase 1 */
  int prOpt_P3_2_3_4_5;     /*equals 1 if the incumbent sol. is proved to be opt. after phase 5                */ 
  int isOpt_P3_2_3_4_5;     /*equals 1 if the incumbent sol. is optimal (not necessarily proven) after phase 5 */
  int prOpt_dom_6;          /*equals 1 if the incumbent sol. is proved to be opt. after phase 6                */ 
  int prOpt_LP_7;           /*equals 1 if the incumbent sol. is proved to be opt. after phase 7                */ 
  int prOpt_Exact_10;       /*equals 1 if the incumbent sol. is proved to be opt. after phase 10               */ 
  int kpHeuristicPhase;     /*number of knapsack problems solved during the heuristic section (phases 2-3) */
  int kpExactPhase;         /*number of knapsack problems solved during the exact section     (phase 10)   */
  int totNodesKpExactPhase; /*total number of nodes (items) in the knapsack problems solved during the exact section */
  int cardP_P3_2_3_4_5;     /*number of nodes in P after phase 5 */
  int cardP_dom_6;          /*number of nodes in P after phase 6 */
  int cardP_LP_7;           /*number of nodes in P after phase 7 */
  int excessSupply;         /*excess supply on the nodes that are used in the optimal solution */
  double opt_val;           /*optimal value obtained */
} Res;

//Contains detailed statistics on the execution of the algorithm
typedef struct GlobalRes{
  uint64_t cpu_sort_0;      /*cpu time spent                   in phase 0  (in microseconds)   */
  uint64_t cpu_ZLP_ZG_1;    /*cpu time spent                   in phase 1  (in microseconds)   */
  uint64_t cpu_P3_2_3_4_5;  /*cpu time spent             in phases 2 to 5  (in microseconds)   */
  uint64_t cpu_dom_6;       /*cpu time spent                   in phase 6  (in microseconds)   */
  uint64_t cpu_LP_7;        /*cpu time spent                   in phase 7  (in microseconds)   */
  uint64_t cpu_filtering_9; /*cpu time spent                   in phase 9  (in microseconds)   */
  uint64_t cpu_exact_10;    /*cpu time spent                   in phase 10 (in microseconds)   */
  uint64_t cpu_tot;         /*total cpu time spent, i.e. in phases 0 to 10 (in microseconds)   */
  int prOpt_ZLP_ZG_1;       /*equals 1 if the incumbent sol. is proved to be opt. after phase 1                */
  int isOpt_ZLP_ZG_1;       /*equals 1 if the incumbent sol. is optimal (not necessarily proven) after phase 1 */
  int prOpt_P3_2_3_4_5;     /*equals 1 if the incumbent sol. is proved to be opt. after phase 5                */ 
  int isOpt_P3_2_3_4_5;     /*equals 1 if the incumbent sol. is optimal (not necessarily proven) after phase 5 */
  int prOpt_dom_6;          /*equals 1 if the incumbent sol. is proved to be opt. after phase 6                */ 
  int prOpt_LP_7;           /*equals 1 if the incumbent sol. is proved to be opt. after phase 7                */ 
  int prOpt_Exact_10;       /*equals 1 if the incumbent sol. is proved to be opt. after phase 10               */ 
  int kpHeuristicPhase;     /*number of knapsack problems solved during the heuristic section (phases 2-3) */
  int kpExactPhase;         /*number of knapsack problems solved during the exact section     (phase 10)   */
  int totNodesKpExactPhase; /*total number of nodes (items) in the knapsack problems solved during the exact section */
  int cardP_P3_2_3_4_5;     /*number of nodes in P after phase 5 */
  int cardP_dom_6;          /*number of nodes in P after phase 6 */
  int cardP_LP_7;           /*number of nodes in P after phase 7 */
  int excessSupply;         /*excess supply on the nodes that are used in the optimal solution           */
  int numberPartials;       /*number of instances for which the optimal solution contains a partial node */
  int nRep;                 /*number of SSFCTPs solved */
  double sum_opt_val;       /*sum of the optimal values of the SSFCTPs solved */
} GlobalRes;

struct GlobalRes* createGlobalRes();
void updateGlobalRes(GlobalRes *gRes, Res * res);
void printGlobalRes(GlobalRes *gRes, int computeSol);

struct Res* create_res();
void cleanMemory(Res* res, P3struct* P3L, Prob* pr, Node* nodes);

void advanceSeed(int * seed);
int initSeed(int seed);
#endif

