#include "combo.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <float.h>
#include <limits.h>

//=============================================== Problem transformation functions ===============================================

// return -1 if e_1 < e_2 and +1 otherwise
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

// sort nodes in non-decreasing order of linearized cost e_j
static void sort_nodes(int n, Node *nodes){
	double max_e = 0.0;
	int i = 0;
	while (max_e <= nodes[i].e)
	{
		max_e = nodes[i++].e;
	}

	if(i < n){ // only call qsort if the nodes were not already sorted
		qsort(nodes, n, sizeof(Node), node_comparator); //Sorting the nodes in non-decreasing order of e_j
	}
	
	for (int i = 0; i < n; ++i){
    	nodes[i].i = i;                           //Node index fits the sorting order
    }                                
            
}

// Transformation of the problem to set all the unit costs c_j to postive value, with min{c_j'} = 0.
// Each unit cost is reduces by theta=min{c_j}. The optimal solution to the modified problem is also the optimal
// solution t the original problem. D*theta must be added to the modified problem's optimal value to obtain the 
// orignal problem's optimal value. This funciton must be called ofter the nodes and problem generation.
// Return D*theta
double setMinCjToZero(Prob *pr, Node *nodes){
	int n = pr->n;
	int D = pr->D;

	double theta = DBL_MAX;
	for (int i = 0; i < n; ++i)
	{
		theta = MIN(theta, nodes[i].c);
	}

	for (int i = 0; i < n; ++i)
	{
		nodes[i].c -= theta;
		nodes[i].e -= theta;
		nodes[i].u -= theta*nodes[i].b;
	}

	pr->lc   = 0;
	pr->uc   -= theta;
	pr->Cmin = pr->lc;
	pr->Cmax = pr->uc;
	pr->LB   = 0; //MAX(DBL_MIN, pr->LB - D*theta);
	pr->UB   = DBL_MAX; //MIN(DBL_MAX, pr->UB - D*theta);

	sort_nodes(n, nodes);
	return D*theta;
}

//============================== Basic functions to produce the data structures that are required by the kaa =============================
// Produce an instance of the data structure "Node" for the desciption of the SSFCTP given in argument
Node* get_nodes(int * b, double * c, double * f, int n)
{
  Node *nodes = (Node*)calloc(n, sizeof(Node));

  int i;
  for(i=0; i<n; i++){
  	nodes[i].i_init = i;
    nodes[i].i = i;
    nodes[i].b = b[i];
    nodes[i].c = c[i];
    nodes[i].f = f[i];
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

  return nodes;
}


// Produce an instance of the data structure "Prob" for the desciption of the SSFCTP given in argument
Prob* get_prob(int D, int * b, double * c, double * f, int n){
  Prob *pr = (Prob*)calloc(1, sizeof(Prob));
  int lb = INT_MAX;
  int ub = INT_MIN;
  double lc = DBL_MAX;
  double uc = DBL_MIN;
  double lf = DBL_MAX;
  double uf = DBL_MIN;

  int i;
  for(i=0; i<n; i++){
    if(b[i] < lb){
      lb = b[i];
    }
    if(c[i] > uc){
      uc = c[i];
    }

    if(c[i] < lc){
      lc = c[i];
    }
    if(c[i] > uc){
      uc = c[i];
    }

    if(f[i] < lf){
      lf = f[i];
    }
    if(f[i] > uf){
      uf = f[i];
    }
  }
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

//================================================= Comparators ==========================================================

// return 1 if lambda_1 > lambda_2, where lambda is the value of lambda used to solve P3(lambda)
// n1 and n2 are pointers to nodes
static int P3struct_comparator(const void *n1, const void *n2)
{
	const P3struct *p1 = (P3struct *)n1;
	const P3struct *p2 = (P3struct *)n2;
	if (p1->lambda < p2->lambda)
		return -1;
	else if (p1->lambda > p2->lambda)
		return +1;
	else
		return 0;
}

// return i iff i<j according to the total order. If i<=j, return 0.
static int i_prec_j(int i, int j, Node * nodes){
	int boolVal = ( (nodes[i].c<nodes[j].c) || (nodes[i].c==nodes[j].c &&  nodes[i].b<nodes[j].b) || (nodes[i].c==nodes[j].c &&  nodes[i].b==nodes[j].b && i<j) );
	return boolVal;
}

// return -1 if n1<n2 according to the article's total order (prop. 11.2), i.e.
// (c_1<c_2) or ((c1=c2) and (b1<b2)) or ((c1=c2) and (b1=b2) and (i1<i2)). Return +1 otherwise. 
// If index_comparator(*n1, *n2)>0, node 1 is excluded from any solution using node 2 as a partial.
// n1 and n2 are pointers to pointers to nodes
static int index_comparator(const void *n1, const void *n2)
{
	const Node *p1 = *(Node **)n1;
	const Node *p2 = *(Node **)n2;
	if ((p1->c<p2->c)||((p1->c==p2->c)&&(p1->b<p2->b))||((p1->c==p2->c)&&(p1->b==p2->b)&&(p1->i<p2->i)))  //n1<n2. i.e. 1 partial => 2 unused
		return -1;
	return +1;
}

// Same comparison rules as in index_comparator in node.c 
// This function in used to insert a single new node in ptrLRnodesC already sorted
static void insertPtrSortedCj(Node ** ptrLRnodesC, int nbNodes, Node newNode) 
{ 
	int i;
	for (i = nbNodes - 1; (i >= 0 && (ptrLRnodesC[i]->c > newNode.c || (ptrLRnodesC[i]->c == newNode.c && ptrLRnodesC[i]->b > newNode.b))); i--) 
		ptrLRnodesC[i + 1] = ptrLRnodesC[i]; 

	ptrLRnodesC[i + 1] = &newNode; 
} 

//============================================ Solve knapsack for fixed partial =======================================================

// Exact resolution of P3|p is partial. 
// This implementation is useful when the partial node p has a unit cost c_p which leads to 
// a negative profit u_j-b_j*c_p in the previous implementation's transformation for node(s) j.
//
// Here, - each node j!=p has weight b_j and profit u_j in the knapsack transformation
//       - p is splitted in ceil(log_2(b_j-LBxp)) nodes with weights corresponding to the powers of two that sum to b_j-LBxp. Each of these nodes i of weight w_i (power of 2) has profit w_i*c_p
//
// The knapsack capacity is given by sum_bj-(D-LBxp)
// Each node j!=p is completely used in the SSFCTP if it in not used in the knapsack and unused in the SSFCTP if it is unused in the knapsack
// LBxp units are used on p, plus the sum of the weights w_i of the "power-of-two nodes" builded after p that are not used in the knapsack.
// 
// Zssfctp = sum_uj - Zkp.
//
// This implementation in more robust than the previous one, since it does not not require the condition (u_j>b_jc_p \forall j!=p) required by the previous implementation
// ** Improved implementation **
//  -> only includes nodes j<p, according to the article's order in the KP
//  -> only includes nodes j with nodes[j].poss_1 in the KP
//  -> directly discards nodes j!=p with nodes[j].poss_1=FALSE in the SSFCTP solution. (nodes that can't be completely used and that are not partial, since p is fixed as a partial, hence they must be unused)
//  -> directly uses nodes j!=p with nodes[j].poss_0=FALSE in the SSFCTP solution. (nodes that can't be unused and that are not partial, since p is fixed as a partial, hence they must be completely used)
//  -> Better decomposition of the remaining offer UBxp-LBxp on the partial node using Martello and Toth decomposition. Now exactly the range {LBxp,LBxp+1,...,UBxp} can be reached by xp
double solveKnapsackFixedPartial(int p, Node * nodes, Prob *pr, int * z, int verbose, int *totNodesKpExactPhase){
	int int_adjustment = 1;
	int n = pr->n;
	double D = pr->D;
	int offer_p = nodes[p].UBx-nodes[p].LBx;
	int nb_nodes_represent_p = 0;
	int weightLast = 0;
	if(offer_p > 0){                                        // in this case, offer_p = nodes[p].UBx-nodes[p].LBx must be decomposed in powers of two
		nb_nodes_represent_p = ceil_log2(offer_p+1);  
		weightLast = offer_p - (pow(2,nb_nodes_represent_p-1)-1);  // the last node of the decomposition has offer b*j = weightLast. Hence, sum_{j=n+1}^{n+d^p} = UBxp-LBxp
	}

	double w[n + nb_nodes_represent_p];          //weight w[i]=b_i
	double profit[n + nb_nodes_represent_p];
	int a[n + nb_nodes_represent_p];
	int idx_ssfctp[n];      //  z[idx_ssfctp[j]] = 1-a[j]. Used to link a KP node to its corresponding SSFCTP node

	int sum_bj_used         = 0; //sum_{j=1}^n(b_j * I[node j is already completely used before the knapsack resolution])
	double sum_uj_used      = 0; //sum_{j=1}^n(u_j * I[node j is already completely used before the knapsack resolution])
	int sum_weight_kp       = 0; //sum_{j=1}^n(b_j * I[node j is included in the knapsack])
	double sum_profit_kp    = 0; //sum_{j=1}^n(u_j * I[node j is included in the knapsack])
	int curIdx              = 0; //index of the knapsack nodes' array where the next node will be added

	for (int i = 0; i < n; ++i) 
	{
		z[i] = 0; //node i is initially unsused in the SSFCTP solution. z[i]=1 if i is completely used and z[i]=0 otherwise (including if i=p)
		if(i==p){
			continue;
		}
		else if(nodes[i].poss_1==0 || partial_node_comparator(&nodes[i], &nodes[p]) > 0){ // i can't be completely used and is not parital or p partial => i unsused
			continue;
		}
		else if(nodes[i].poss_0==0){ // i can't be unused and is not parital => i completely used
			z[i] = 1;
			sum_bj_used += nodes[i].b;
			sum_uj_used += nodes[i].u;
		}
		else{ // i can be unused or completely used. This will be decided when solving the knapsack
			idx_ssfctp[curIdx] = i; //the node stocked in position curIdx to solve the KP is node i. 
			a[curIdx] = 0;
			profit[curIdx] = nodes[i].u;
			w[curIdx] = (double)nodes[i].b;
			sum_weight_kp += nodes[i].b;
			sum_profit_kp += nodes[i].u;
			curIdx += 1;
		}	
	}

	int idx_decomposed_p = curIdx; //next index to add a node from p's decomposition 
	if(nb_nodes_represent_p>0){
		for (int i = 0; i < nb_nodes_represent_p-1; ++i) // division of UBxp-LBxp in powers of two 
		{
			int weight = pow(2, i); //ex: if offer_p = 19, nb_nodes_represent_p = 5. The first weight produced is 2^(5-0-1) = 16
			a[idx_decomposed_p] = 0;
			profit[idx_decomposed_p] = weight*nodes[p].c;
			w[idx_decomposed_p] = (double)weight;
			sum_weight_kp += weight;
			sum_profit_kp += profit[idx_decomposed_p];
			idx_decomposed_p += 1;
		}
		a[idx_decomposed_p] = 0;
		profit[idx_decomposed_p] = weightLast*nodes[p].c;
		w[idx_decomposed_p] = (double)weightLast;
		sum_weight_kp += weightLast;
		sum_profit_kp += profit[idx_decomposed_p];
		idx_decomposed_p += 1;
	}
	

	int64_t C = (int64_t)(sum_weight_kp - (D-(sum_bj_used + nodes[p].LBx)));
	if(verbose==2){
		printf("node %d partial, number of nodes passed to combo: %d\n", p, idx_decomposed_p);
	}
	*totNodesKpExactPhase += idx_decomposed_p;

	if(C < 0 || sum_bj_used + nodes[p].LBx > D){ //no feasible solution with p partial respects the already known (valid) constraints on the nodes' usage
			printf("*********************************************\n");
			printf("* NO FEASIBLE SOLUTION WITH NODE %d PARTIAL *\n",p); 
			printf("*********************************************\n");
		return pr->UB+1;
	}

	double Zssfctp; //objective function of the corresponding SSFCTP solution
	double M = sum_uj_used + nodes[p].f + nodes[p].LBx*nodes[p].c + sum_profit_kp; //  sum_uj_used + nodes[p].f + (nodes[p].LBx - extra)*nodes[p].c is already paid. sum_profit_kp-Zkp is added, where Zkp is the sum of cost u_j of nodes putted in the optimal knapsack
	double lb_combo = (M - pr->UB);
	double ub_combo = (M - nodes[p].LB);

	if(sum_weight_kp <= C){ // if all the offer fits in the KP, the solution is trivial. Everything is placed in the kp. 
		for (int i = 0; i < idx_decomposed_p; ++i)
		{
			a[i]=1;
		}
	}
	else{ // non-trivial solution
		combo_solve(profit, w, a, C, lb_combo, ub_combo, idx_decomposed_p);
	}

	Zssfctp = sum_uj_used;

	for (int i = 0; i < curIdx; ++i)
	{
		z[idx_ssfctp[i]] = 1-a[i];     // Node idx_ssfctp[i] is completely use if a[i]=0 (has not been taken in the KP) and is unused if a[i]=1 (has been taken in the KP)
		Zssfctp += (1-a[i])*profit[i]; // (1-a[i])*profit[i] == z[idx_ssfctp[i]]*nodes[idx_ssfctp[i]].u
	}

	int units_p_unused = 0; //number of units of p taken in the KP (unused in the SSFCTP). If no unit is removed (units_p_unused = 0), hence xp=UBx
	for (int i = curIdx; i < idx_decomposed_p; ++i)
	{
		units_p_unused += a[i]*w[i];
	}
	int xp = nodes[p].UBx-units_p_unused;
	Zssfctp += nodes[p].f + xp*nodes[p].c;

	if(Zssfctp > pr->UB){
		Zssfctp = pr->UB+1;
		if(verbose==2){
			printf("P3|%d is partial > UB\n\n", p); //Here we do not explicitly return the value of Zssfctp calculated, because it does not necessarily represent the real value of P3|p is partial, 
		}                                           //since when (P3|p is partial)>UB, combo returns all the values a[i]=0, hence no valid solution. All we have in this situation is the proof that (P3|p is partial)>UB.
	}                                              
	else{
		if(verbose==2){
			printf("P3|%d is partial = %f\n\n", p, Zssfctp);
		}
	}
	return Zssfctp;
}


// Exact resolution of P3 using a call of combo for each potential partial node
void solveKnapsackFixedPartial_all(Node * nodes, Prob *pr, int *totCallsComboExactPhase, int verbose, int *totNodesKpExactPhase, int computeSol){
	if(pr->nb_poss_p == 0)return;

	for (int p = 0; p < pr->n; ++p)
	{
		if(nodes[p].poss_p){
			if(nodes[p].LB >= pr->UB){
				nodes[p].poss_p=0;
				pr->nb_poss_p-=1;
			}
			else{
				*totCallsComboExactPhase +=1;
				int z[pr->n];
				double Zssfctp = solveKnapsackFixedPartial(p, nodes, pr, z, verbose, totNodesKpExactPhase);
				if(Zssfctp < pr->UB){
					pr->UB = Zssfctp;

					//=== Update the incumbent solution ===
					if(computeSol){
						int offer_comp = 0;
						for (int i = 0; i < pr->n; ++i){ 
							int x_i = z[i]*nodes[i].b;
							offer_comp += x_i;
		    				pr->incumbent[i] = x_i; 
		    			}
		    			pr->incumbent[p] = pr->D - offer_comp;
		    		}
	    			//======================================

				}
				nodes[p].poss_p=0;
				pr->nb_poss_p-=1;
			}
		}
	}
	pr->LB = pr->UB;
}

//================================================= 10-10.1: P3 calculations ==========================================================
//Use known values of P3(lambda) to get an upper bound on P3(lambda') for a new value lambda' //lb_ub must be an array of 2 doubles containing -1 in each case
int bounds_P3lambda(double* lb_ub, double lambda, P3struct * P3L, int sizeP3L){//result in lb_ub, return the incumbent solution's index from P3L for combo(lambda)
	int index_incumbent = -1; //no incumbent solution yet
	if(P3L[0].P3lambda == 0){//L in empty. No bound can be calculated
		return index_incumbent;
	}
	if(lambda <= P3L[0].lambda){//lambda < P3L[0].lambda
		lb_ub[0] = P3L[0].P3lambda;
		lb_ub[1] = P3L[0].P3lambda - P3L[0].surplus*(lambda - P3L[0].lambda);
		index_incumbent = 0;
		return index_incumbent;
	}
	else{
		for (int i = 1; i < sizeP3L; ++i)
		{
			if(P3L[i].lambda>=lambda){//P3L[i-1].lambda < lambda < P3L[i].lambda
				lb_ub[0] = P3L[i-1].lambda - (lambda - P3L[i-1].lambda)*(P3L[i-1].P3lambda-P3L[i].P3lambda)/(P3L[i].lambda-P3L[i-1].lambda); //line pasing through points (lambda_1, P3(lambda_1)), (lambda_2, P3(lambda_2)), evaluated in x=lambda
				double ub = P3L[i-1].P3lambda - P3L[i-1].surplus*(lambda - P3L[i-1].lambda); //UB given by the opt. solution to P3(lambda_1) 
				index_incumbent = i-1;
				double ub2 = P3L[i].P3lambda - P3L[i].surplus*(lambda - P3L[i].lambda);      //UB given by the opt. solution to P3(lambda_2) 
				if(ub2 < ub){
					index_incumbent =i;
					ub = ub2;
				}
				lb_ub[1]=ub;
				return index_incumbent;
			}
			if(P3L[i].P3lambda == 0){ //no lambda in L bigger than L (and less than sizeP3L values in L) (no lower bound => lb=0)
				lb_ub[0] = 0;
				lb_ub[1] = P3L[i-1].P3lambda - P3L[i-1].surplus*(lambda - P3L[i-1].lambda); 
				index_incumbent = i-1;
				return index_incumbent; //no incumbent solution
			}
		}
	}
	//no lambda in L bigger than L (and exactly sizeP3L values in L) (no lower bound -> lb=0)
	lb_ub[0] = 0;
	lb_ub[1] = P3L[sizeP3L].P3lambda - P3L[sizeP3L].surplus*(lambda - P3L[sizeP3L].lambda); 
}

// returns a "P3struct" data structure with lambda, cp, p, P3(lambda)
// Zssfctp(z*) and surplus, with z* the optimal solution to P3(lambda).
// the array z passed in argument is modified to z* after the execution
// of the function
P3struct * P3(double lambda, Node * nodes, int n, int D, int * z, P3struct * P3L, int sizeP3L){
	int int_adjustment = 1;
	P3struct * p3 = (P3struct*)calloc(1, sizeof(P3struct));
	double w[n];
	double profit[n];  //profit profit[i]=u_i-b_i*lambda on node i
	int a[n];

	int sumB = 0;
	double sumProfit = 0;

	for (int i = 0; i < n; ++i)
	{
		sumB += nodes[i].b;
		a[i] = 0;
		profit[i] = nodes[i].u - lambda*nodes[i].b;
		sumProfit += profit[i];
		w[i] = (double)nodes[i].b;
	}

	int64_t C = (int64_t)(sumB - D);
	
	double lb_combo, ub_combo;
	double lb_ub_P3lambda[2] = {-1,-1};
	int index_incumbent = bounds_P3lambda(lb_ub_P3lambda, lambda, P3L, sizeP3L);
	if(lb_ub_P3lambda[1] < 0){ //no ub (nor lb) could be calculated for P3(lambda), since L is empty
		lb_combo = 0;
		ub_combo = lambda*D+sumProfit;
	}
	else{
		lb_combo = lambda*D + sumProfit - lb_ub_P3lambda[1]; //lb_ub_P3lambda[1] = ub(P3(lambda))
		ub_combo = lambda*D + sumProfit - lb_ub_P3lambda[0]; //lb_ub_P3lambda[0] = lb(P3(lambda))
	}

	double combolambda;
	double P3lambda;

	if(lb_combo < ub_combo){                                                    
		double Z_combo = combo_solve(profit, w, a, C, floor(lb_combo), ceil(ub_combo), n); //TODO, remove floor() and ceil(), (still required because of a Segmentation fault)
		if(Z_combo <= lb_combo){ //no feasible solution to combo with objective value better than lb_combo exists
			p3->lambda = lambda;
			p3->cp = P3L[index_incumbent].cp;
			p3->p = P3L[ index_incumbent].p;
			p3->surplus = P3L[ index_incumbent].surplus;
			p3->Zssfctp = P3L[ index_incumbent].Zssfctp;
			p3->P3lambda = lb_ub_P3lambda[1];
			return p3;
		}
		else{
			combolambda = 0;
			for (int i = 0; i < n; ++i)
			{
				combolambda+=a[i]*profit[i]; //use the double value (full precision)
			}
		}
	}
	else{
		p3->lambda = lambda;
		p3->cp = P3L[index_incumbent].cp;
		p3->p = P3L[ index_incumbent].p;
		p3->surplus = P3L[ index_incumbent].surplus;
		p3->Zssfctp = P3L[ index_incumbent].Zssfctp;
		p3->P3lambda = lb_ub_P3lambda[1];
		return p3;
	}
	P3lambda = lambda * D + sumProfit - combolambda;

	double cp = 0;
	int p = 0;
	int offer = 0;

	for (int i = 0; i < n; ++i)
	{
		z[i] = 1-a[i];
		if(z[i]*nodes[i].c > cp){
			cp = nodes[i].c;
			p = i;
		}
		offer += z[i]*nodes[i].b;
	}
	int surplus = offer - D;
	p3->surplus  = surplus;
	if(surplus <= nodes[p].b){
		p3->Zssfctp = P3lambda + (lambda-cp)*surplus; //on remet le remboursement de lambda*surplus et on recoit un remboursement de cp*surplus à la place
	}
	else{
		p3->Zssfctp = P3lambda + (lambda)*surplus - nodes[p].u ; //on retire p de la solution complètement et on remet le remboursement pour le surplus
	}                                                            //ceci donne une UB, car on reste avec une solution où on paye des unités pour rien

	p3->lambda = lambda;
	p3->cp = cp;
	p3->p = p;
	p3->P3lambda = P3lambda;
	return p3;
}


// returns an array of "P3struct" data structures for each lambda in L
// the reseach is implemented as explained in the article
// Cmin can be imporved during the execution
// Use max_{j in {1,2,...,s}}{c_j} as initial value of lambda in P3Search if true, and min{e_s, max_{j in {1,2,...,n}}{c_j}} otherwise
P3struct * P3search(Prob *pr, Node *nodes,  int sizeP3L, int *totCallsComboP3Phase, int computeSol, int alternative_lambda_init){
    //First lambda is min{e_s, max_{i<=n}(c_i)}. In the most frequent case (i.e. e_s >= max{c_j}),
    //When e_s < max{c_j}, to avoid uninformative calls to
    //combo (which becomes a max-kp with less than C units of weight on nodes with positive profit, which has a trivial solution
    //that consist in taking all the positive profit nodes and not taking the others/ this may lead to "backtracking problems" in combo),
    //instead of calling P3(max_cj), we first call P3(e_s) and then perform the search as usual, except that we stop the search whenever 
    //lambda'>es, i.e. we try to do a search on the forbidden space. 
    //In the latter case (es<max{cj}), P3 does not lead to a global lower bound on the objective value, hence does not allow to conlcude 
    //that the optimal has been found directly. However, it still gives a valid upper bound (generally optimal) and a lower bound on
    //[P3|p partial] for each node p with c_p<=e_s.
	int offer = 0;
	double es = 0;
	double max_cj = 0;
	double lambda = 0.0; 
	double max_cj_j_in_1_s = 0; 
	for (int i = 0; i < pr->n; ++i)
	{
		if(offer < pr->D){
			offer += nodes[i].b;
			if(offer >= pr->D){
				es = nodes[i].e;
			}
			if(nodes[i].c > max_cj_j_in_1_s){
				max_cj_j_in_1_s = nodes[i].c;
			}
		}
		max_cj = MAX(max_cj, nodes[i].c);
	}
	double lambda_max = MIN(max_cj, es);

	if(alternative_lambda_init){
		lambda = max_cj_j_in_1_s; //alternative initial lambda:  max_{j in {1,2,...,s}}{c_j} (aims to identify the optimal solution quickly, does not lead to a global LB)
	}
	else{
		lambda = lambda_max; //the first value of lambda used is lambda_max:=min{e_s, max_{i<=n}(c_i)} (leads to a global LB if max_{i<=n}(c_i) <= e_s, but is less likely to identify the optimal solution directly)
	}                                                             

	P3struct * P3L = (P3struct*)calloc(sizeP3L, sizeof(P3struct)); // maximum sizeP3L calls to P3(lambda)
	double L[sizeP3L];
	int idxlambda = 0;
	P3struct * p3 = (P3struct*)calloc(1, sizeof(P3struct)); // current P3(lambda)
	int lambda_in_L = 0;
	int z[pr->n];

	double max_in_L; // maximum value of lambda for which P3(lambda) has been computed

	while(lambda_in_L == 0 && idxlambda <= sizeP3L-2 && lambda<=lambda_max){ // a maximum of sizeP3L-1 values of lambda are calculated here
			//printf("lambda = %lf\n", lambda); //TODO REMOVE AFTER FINISHING TESTING FOR INITIAL VALUE OF LAMBDA
			max_in_L = MAX(lambda, max_in_L);
	    p3 = P3(lambda, nodes, pr->n, pr->D, z, P3L, sizeP3L);       // new value of P3(lambda)
	    if(p3->Zssfctp <= pr->UB){
	    	pr->UB = p3->Zssfctp;

	    	//=== Update the incumbent solution ===
	    	if(computeSol){
				int offer_comp = 0;
				int p = -1;
				for (int i = 0; i < pr->n; ++i){ 
					int x_i = z[i]*nodes[i].b;
					offer_comp += x_i;
					pr->incumbent[i] = x_i; 
					if(x_i > 0 && (p==-1 || nodes[i].c > nodes[p].c || (nodes[i].c >= nodes[p].c && nodes[i].b >= nodes[p].b) )){
						p = i;
					}
				}
				pr->incumbent[p] -= (offer_comp - pr->D);
			}
			//======================================
	    }
	    P3L[idxlambda] = *p3;
	    L[idxlambda] = p3 -> lambda;
	    if(p3 -> lambda == p3 -> cp){
	    	if(pr->Cmin < p3 -> lambda){ //P3(lambda) = Zssfctp(z*), where z* is the optimal solution to P3(lambda).
	    		pr->Cmin = p3 -> lambda; //Hence any solution to P3 without partial or with partial p with cp<=lambda has objective value >= Zssfctp(z*)
	    		pr->cmin_closed_interval = 0; //we then only search for partials with c_j>lambda. (open interval)
	    	}
	    }

	    idxlambda++;
	    lambda = p3 -> cp;

	    for (int i = 0; i < idxlambda; ++i)
	    {
	    	if(lambda == L[i]){
	    		lambda_in_L = 1;
	    	}
	    }

	    // also stop the search if the optimal solution to P3(lambda) does not include a partial node (i.e. the offer on the nodes that are used in the solution is exactly equal to D)
	    if(p3->surplus == 0)
	    	break;
	}

	qsort(P3L, idxlambda, sizeof(P3struct), P3struct_comparator);  //Sorting the P3struct in increasing order of lambda

	if(pr->UB > P3L[0].P3lambda){  // if there is no value lambda \in L such that Z_{UB} <= Z^{KP}_{\lambda}, we solve P3(0) so we will find the optimal solution if it does not include a partial node
		printf("P3(0) has to be solved after the execution of algorithm P3Search\n");
		p3 = P3(0, nodes, pr->n, pr->D, z, P3L, sizeP3L);      
	    if(p3->Zssfctp <= pr->UB){
	    	pr->UB = p3->Zssfctp;
	    	//=== Update the incumbent solution ===
	    	if(computeSol){
				int offer_comp = 0;
				int p = -1;
				for (int i = 0; i < pr->n; ++i){ 
					int x_i = z[i]*nodes[i].b;
					offer_comp += x_i;
					pr->incumbent[i] = x_i; 
					if(x_i > 0 && (p==-1 || nodes[i].c > nodes[p].c || (nodes[i].c >= nodes[p].c && nodes[i].b >= nodes[p].b) )){
						p = i;
					}
				}
				pr->incumbent[p] -= (offer_comp - pr->D);
			}
			//======================================
	    }
	    P3L[idxlambda] = *p3;
	    idxlambda++;
	    qsort(P3L, idxlambda, sizeof(P3struct), P3struct_comparator);  //Sorting the P3struct in increasing order of lambda (to insert P3(0) at the beginning)
	}

	if(max_in_L == max_cj){ //P3(max_cj) has been calculated. P3(max_cj) is a valid LB on the objective value and LB(P3|p partial) can be calulated for each node p 

		if(P3L[idxlambda-1].P3lambda == P3L[idxlambda-1].Zssfctp){ //P3(lambda) is a lower bound for any solution with partial p with cp<=lambda.
			  pr->LB = pr->UB;                                       //Here lambda==max{c_j} => P3(lambda) is a lower bound for any solution with a partial node
		    pr->Cmin = P3L[idxlambda-1].lambda;                    //Also, the optimal solution to P3 is a feasible solution to P3, with the same objective value (P3(lambda)=Zssfctp(z*)), where z* is the optimal solution to P3(lambda) 
		    pr->cmin_closed_interval = 0;                          //same justification as above
		    pr->Cmax = P3L[idxlambda-1].lambda;
		    pr->cmax_closed_interval = 0;
		}
		else{
	    	pr->LB = MAX(P3L[idxlambda-1].P3lambda,pr->LB); // Global lower bound given by P3(max{c_j})                             
		}
	}

	//else{ //(max_lambda = e_s < max{c_j}) special case: No global lower bound can be calculated. 
		//No global lower bound can be calculated. 
	//}

	//Update on Cmin if UB_lambda_eq_cp = FALSE (Explained in article sec. 3.2.2)
	int UB_lambda_eq_cp = 0;  // TRUE if the current UB has been obtained with a P3(lambda) : lambda==cp (Cmin directly given by cp)
	for (int i = 0; i < idxlambda; ++i){
		if(P3L[i].lambda == P3L[i].cp && P3L[i].Zssfctp <= pr->UB){  //*(1+EPSILON)) {
			pr->Cmin = MAX(pr->Cmin, P3L[i].cp);
			UB_lambda_eq_cp = 1;
			break;
		}
	}
	if(UB_lambda_eq_cp == 0){
		int idx0;
		int idx1 = 0;
		for (int i = 1; i < idxlambda; ++i)
		{
			if(P3L[i].P3lambda < pr->UB*(1-EPSILON)){
				idx1 = i;
				idx0 = i-1;
				break;
			}
		}
		if(idx1 > 0){ // means that there are two multipliers lambda0 and lambda1 allowing to compute the bound Cmin
			double lambda0   = P3L[idx0].lambda; 
			double P3lambda0 = P3L[idx0].P3lambda; 
			double lambda1   = P3L[idx1].lambda; 
			double P3lambda1 = P3L[idx1].P3lambda; 
			pr->Cmin = lambda0 + (lambda1-lambda0)*(P3lambda0-pr->UB)/(P3lambda0-P3lambda1);
			pr->cmin_closed_interval = 0; //the incumbent can't be strictly improved with a partial having a unit cost <=cmin
		}
	}

	for (int i = 0; i < pr->n; ++i) //Update the list of possibly partial nodes.
	{
		if(nodes[i].c<pr->Cmin ||(nodes[i].c==pr->Cmin && pr->cmin_closed_interval==0)){
			nodes[i].LB = pr->UB+1;
			nodes[i].poss_p = 0;
			pr->nb_poss_p-=1;
		}
	}

	*totCallsComboP3Phase += idxlambda;
	return P3L;
}

//================================================= 10.2: LB on (P3|p is partial) ==========================================================
//LBp1: implementation in O(|L|) for each node (linear search)
double LBp1(int p, Node * nodes, P3struct * P3L, int sizeP3L){
	double cp = nodes[p].c;
	for (int i = 0; i < sizeP3L; ++i)
	{
		if(cp == P3L[i].lambda) return P3L[i].P3lambda;
		if(cp < P3L[i].lambda){
			double lambda0   = P3L[i-1].lambda; 
			double P3lambda0 = P3L[i-1].P3lambda; 
			double lambda1   = P3L[i].lambda; 
			double P3lambda1 = P3L[i].P3lambda; 
			return P3lambda0 + (P3lambda1-P3lambda0)/(lambda1-lambda0)*(cp-lambda0);
		}
	}
}

// Calculate LBp1 on each potential partial node
void LBp1_all(Node * nodes, P3struct * P3L, Prob * pr, int n, int sizeP3L){
	if(pr->nb_poss_p == 0)return;

	double minLBp = pr->UB;
	for (int i = 0; i < n; ++i)
	{
		if(!nodes[i].poss_p) continue;
		nodes[i].LB = LBp1(i, nodes, P3L, sizeP3L);
		if(nodes[i].LB < minLBp) minLBp = nodes[i].LB;
	}
	pr->LB = MAX(minLBp, pr->LB);
}

//Solve the unconstrained LR. Returns the corresponding LRstruct. After the execution, LRnodesC contains the 
//indexes of the nodes {0,1,2,...,s} used in the linear relaxation, sorted by non-decreasing order of c_j
LRstruct * baseLR(Node * nodes, Prob * pr, int n, Node ** ptrLRnodesC){
	LRstruct * LR = (LRstruct*)calloc(1, sizeof(LRstruct));
	int Dem = pr->D;

	for (int i = 0; i < n; ++i)
	{
		if(nodes[i].b < Dem){
			LR->cost += nodes[i].u;
			Dem -= nodes[i].b;
		}
		else{
			LR->cost += Dem*nodes[i].e;
			LR->s = i;
			LR->xs = Dem;
			Dem = 0;
			break;
		}
	}
	qsort(ptrLRnodesC, LR->s+1, sizeof(ptrLRnodesC[0]), index_comparator);
	return LR;
}

//remove currentNodeLL from the linkedList. currentNodeLL is now currentNodeLL->prev
static void removeFromLL_backwards(NodeLinkedList **currentNodeLL){
	if((*currentNodeLL)->prev == NULL){
		(*currentNodeLL)->next->prev=NULL; //no more nodes available
		(*currentNodeLL) = NULL;
	}
	else{                               //general case. currentNodeLL removed from the linked list
		NodeLinkedList *newCurrent = (*currentNodeLL)->prev;
		newCurrent->next = (*currentNodeLL)->next;
		(*currentNodeLL)->next->prev = newCurrent;
		//free(currentNodeLL);//TODO ENLEVER SI BUG
		(*currentNodeLL) = newCurrent;
	}
}

//remove currentNodeLL from the linkedList. currentNodeLL is now currentNodeLL->next
static void removeFromLL_forwards(NodeLinkedList **currentNodeLL){
	if((*currentNodeLL)->next == NULL){
		(*currentNodeLL)->prev->next=NULL; //no more nodes available
		(*currentNodeLL) = NULL;
	}
	else{                               //general case. currentNodeLL removed from the linked list
		NodeLinkedList *newCurrent = (*currentNodeLL)->next;
		newCurrent->prev = (*currentNodeLL)->prev;
		(*currentNodeLL)->prev->next = newCurrent;
		//free(currentNodeLL);//TODO ENLEVER SI BUG
		(*currentNodeLL) = newCurrent;
	}
}

// Calculate LBp2(strong relaxation), LBxp2(strong relaxation) and UBxp2(strong relaxation) on each potential partial node: : O(n(logn + min(n, max_{b_j}/min_{b_j}))) implementation
void LBp2_LB_UBxp2_all_efficient(Node * nodes, Prob * pr, int n){
	if(pr->nb_poss_p == 0)return;

	// nodesLR: Max-heap (using total order from the article) containing the nodes currently used in the LR 
	Node* nodesLR[n];
	int nbElementsLR = 0;
	int s = 0;
	double ZLR = 0;
	int Dem = 0;
	build_heap_LR(pr->n, pr->D, &ZLR, &s, nodes, nodesLR, &nbElementsLR); //O(n logn)

	// firstNodeLL: Linked list of pointers to nodes that can still be used with the next partial nodes to be tried,
	//              sorted in non-decreasing order of linearized cost c_j
	NodeLinkedList *firstNodeLL, *splitNodeLL, *lastNodeLL;
	build_linkedList_nodes(&firstNodeLL, &splitNodeLL, &lastNodeLL, nodes, s, n); //O(n)
	NodeLinkedList *currentNodeLL = splitNodeLL;

	// nodesPart: Array list of pointers to nodes that may be partial in an optimal solution to P3, sorted in 
	//			  increasing order (using total order from the article)
	Node *nodesPart[n];
	int nb_poss_part;
	build_arraylist_partials(n, nodes, nodesPart, &nb_poss_part); //O(n logn)

	// for each potential partial node p, calculate LBp and (LBxp or UBxp)
	for (int p = nb_poss_part-1; p >= 0; --p){

		// *** Phase 1: remove unusable nodes from the LR ***
		while(partial_node_comparator(nodesPart[p], nodesLR[0]) < 0){ // remove from the LR all nodes that
			Dem += nodesLR[0]->x;                                     // can't be used when p is partial
			ZLR -= nodesLR[0]->x*nodesLR[0]->e; 
			nodesLR[0]->x=0;
			heap_pop(nodesLR, &nbElementsLR);
		}

		// *** Phase 2: insert new nodes in the LR to fill the demand ***
		currentNodeLL = splitNodeLL; //insert new nodes in the LR to fill the demand
		while(Dem > 0){              //starting from node splitNodeLL->ptrNode
			if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){ //currentNodeLL->ptrNode can't be used
				removeFromLL_forwards(&currentNodeLL);
    			if(currentNodeLL == NULL){ //Not enough offer to fill the demand using only nodes usables with p (or any following partial node to be tested)
    				printf("No feasible solution with partials <= nodesPart[p] (<= according to the article's ordering)\n"); 
    				for (int q = p; q >= 0; --q)
    				{
    					nodesPart[q]->LB = pr->UB+1;
    					nodesPart[q]->poss_p = 0;
    					pr->nb_poss_p -= 1;
    				}
    				//free(*nodesLR);   // free the heap
    				//	*nodesLR = NULL;
			    	//free(*nodesPart); // free the arraylist 
			    	//	*nodesPart = NULL;
			    	util_free_LL(firstNodeLL, lastNodeLL); // free the linkedlist
			    	return;
			    }
    			continue; //continue on the next node
    		}
			else{  //currentNodeLL->ptrNode will be used
				if(currentNodeLL->ptrNode->x > 0){  // currentNodeLL->ptrNode was already used in the LR (as the split node)
					int remOfferSplit = currentNodeLL->ptrNode->b - currentNodeLL->ptrNode->x; //remaining offer n the split node
					if(remOfferSplit >= Dem ){ //can fill the demand with the split node
						currentNodeLL->ptrNode->x += Dem;
						ZLR += Dem*currentNodeLL->ptrNode->e;
						Dem = 0;
					}
					else{ //cannot fill the demand with the split node. Use it completely
						currentNodeLL->ptrNode->x = currentNodeLL->ptrNode->b;
						ZLR += remOfferSplit*currentNodeLL->ptrNode->e;
						Dem -= remOfferSplit;
					}
				}
				else{ // currentNodeLL->ptrNode in unused for now. It must be added to the LR
					if(currentNodeLL->ptrNode->b >= Dem ){ //can fill the demand with this node
						currentNodeLL->ptrNode->x = Dem;
						ZLR += Dem*currentNodeLL->ptrNode->e;
						Dem = 0;
					}
					else{ //cannot fill the demand with this node. Use it completely
						currentNodeLL->ptrNode->x = currentNodeLL->ptrNode->b;
						ZLR += currentNodeLL->ptrNode->b*currentNodeLL->ptrNode->e;
						Dem -= currentNodeLL->ptrNode->b;
					}
					heap_insert(nodesLR, currentNodeLL->ptrNode, &nbElementsLR); //insert the new node in the LR heap
				}
			}
			if(Dem > 0){
				currentNodeLL = currentNodeLL->next; //continue on the next node
			}
			else{
				splitNodeLL = currentNodeLL; //new split node (Phase 2 is done)
			}
		}

		// *** Phase 3: adjust the solution after fixing p as a partial node ***
		double LBp2 = ZLR + nodesPart[p]->f - nodesPart[p]->x*(nodesPart[p]->e - nodesPart[p]->c); // pay the real cost for x units on p instead of the linearized cost
		int xp, xs;//__PHASE4 data                                                                 // units used on phase 3's adjusted LR on p and the split node (noted s) of the adjusted solution
		if(nodesPart[p]->c > currentNodeLL->ptrNode->e){     //CASE 1) c_p > e_s.
			if(nodesPart[p]->x < nodesPart[p]->LBx){ // A) Must respect the usage bounds. Adjust the LR by transferring units from s,s-1,.. to p (unprofitable)
				//printf("1-A\n");
				int missingUnitsPartial = nodesPart[p]->LBx - nodesPart[p]->x;   //missingUnitsPartial MUST be added to p. i.e. commanded to p
				for(;;){
					if(currentNodeLL == NULL){ // Impossible to respect LBxp
						LBp2 = pr->UB+1; 
						break;
					}
					if(currentNodeLL->ptrNode->x == 0){ // this means that the current node has been removed from the LR (not usable anymore with the remaining partials to be tried)
						removeFromLL_backwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->prev
					}
                    else if(currentNodeLL->ptrNode->x < missingUnitsPartial){ //transfer all the offer from currentNodeLL to p. Not sufficient yet to respect the bound
                    	LBp2 += currentNodeLL->ptrNode->x * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	missingUnitsPartial -= currentNodeLL->ptrNode->x;
                        currentNodeLL = currentNodeLL->prev;    //will next remove units from a less expensive node (s-1,...)
                    }
                    else{ //add missingUnitsPartial to p to respect LBxp exactly. (Phase 3 is done)
                    	LBp2 += missingUnitsPartial * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	xs = currentNodeLL->ptrNode->x - missingUnitsPartial;//__PHASE4 data
                    	missingUnitsPartial = 0;
                    	xp = nodesPart[p]->LBx;//__PHASE4 data
                    	break;
                    }
                }
            }
			else{                                    // B) General case. Adjust the LR by transferring units from p to s,s+1,... (profitable)
				//printf("1-B\n");
				int remUnitsPartial = nodesPart[p]->x - nodesPart[p]->LBx;       //UP TO remUnitsPartial may be removed from p. i.e. removed from the command to p
				for(;;){
					int remOfferCurrent = currentNodeLL->ptrNode->b - currentNodeLL->ptrNode->x; //up to remOfferCurrent units may be transfered from p to currentNodeLL
					if(nodesPart[p]->i == currentNodeLL->ptrNode->i){ //skip currentNodeLL if it is the same node as p
						if((currentNodeLL->next) == NULL || (currentNodeLL->next)->ptrNode->e >= nodesPart[p]->c){ //No longer profitable (or possible) to move units from p to s+1... (Phase 3 is done)
							xp = nodesPart[p]->LBx + remUnitsPartial;//__PHASE4 data
							xs = xp;//__PHASE4 data
							break;
						}
						currentNodeLL = currentNodeLL->next;    //will next add units on a more expensive node (s+1,...)
					}
					if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){// this means that the current node is not usable anymore with the remaining partials to be tried
						removeFromLL_forwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->next
					}
                    else if(remOfferCurrent < remUnitsPartial){ //transfer remOfferCurrent units from p to currentNodeLL. More units still can be transfered from p to s+1,...
                    	LBp2 -= remOfferCurrent * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	remUnitsPartial -= remOfferCurrent;
                        if((currentNodeLL->next) == NULL || (currentNodeLL->next)->ptrNode->e >= nodesPart[p]->c){ //No longer profitable (or possible) to move units from p to s+1... (Phase 3 is done)
                        	xs = currentNodeLL->ptrNode->b;//__PHASE4 data
							xp = nodesPart[p]->LBx + remUnitsPartial;//__PHASE4 data
							break;
						}
						currentNodeLL = currentNodeLL->next;    //will next add units on a more expensive node (s+1,...)
					}
                    else{ //remove remUnitsPartial from p to respect LBxp exactly (Phase 3 is done)
                    	LBp2 -= remUnitsPartial * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	xs = currentNodeLL->ptrNode->x + remUnitsPartial;//__PHASE4 data
                    	remUnitsPartial = 0;
                    	xp = nodesPart[p]->LBx;//__PHASE4 data
                    	break;
                    }
                }
            }
        }
		else if(nodesPart[p]->c < currentNodeLL->ptrNode->e){//CASE 2) c_p < e_s.
			if(nodesPart[p]->x > nodesPart[p]->UBx){ // A) Must respect the usage bounds. Adjust the LR by transferring units from p to s,s+1,... (unprofitable)
				//printf("2-A\n");
				int exceedingUnitsPartial = nodesPart[p]->x - nodesPart[p]->UBx; //exceedingUnitsPartial MUST be removed from p
				for(;;){
					if(nodesPart[p]->i == currentNodeLL->ptrNode->i){ //skip currentNodeLL if it is the same node as p
						currentNodeLL = currentNodeLL->next;    //will next add units on a more expensive node (s+1,...)
					}
					if(currentNodeLL == NULL){ // Impossible to respect UBxp
						LBp2 = pr->UB+1; 
						break;
					}
					int remOfferCurrent = currentNodeLL->ptrNode->b - currentNodeLL->ptrNode->x; //up to remOfferCurrent units may be transfered from p to currentNodeLL
					if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){// this means that the current node is not usable anymore with the remaining partials to be tried
						removeFromLL_forwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->next
					}
					else if(remOfferCurrent < exceedingUnitsPartial){ //transfer remOfferCurrent units from p to currentNodeLL. More units still MUST be transfered from p to s+1,...
						LBp2 += remOfferCurrent * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
						exceedingUnitsPartial -= remOfferCurrent;
                        currentNodeLL = currentNodeLL->next;    //will next add units on a more expensive node (s+1,...)
                    }
                    else{ //remove exceedingUnitsPartial from p to respect UBxp the bound exactly. (Phase 3 is done)
                    	LBp2 += exceedingUnitsPartial * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
                    	xs = currentNodeLL->ptrNode->x + exceedingUnitsPartial;//__PHASE4 data
                    	exceedingUnitsPartial = 0;
                    	xp = nodesPart[p]->UBx;//__PHASE4 data
                    	break;
                    }
                }
            }
			else{                                    // B) General case. Adjust the LR by transferring units from s,s-1,.. to p (profitable)
				//printf("2-B\n");
				int remOfferPartial = nodesPart[p]->UBx - nodesPart[p]->x;       //UP TO remOfferPartial may be added to p
				for(;;){
					if(nodesPart[p]->i == currentNodeLL->ptrNode->i){ //skip currentNodeLL if it is the same node as p
						if((currentNodeLL->prev) == NULL || (currentNodeLL->prev)->ptrNode->e <= nodesPart[p]->c){ //No longer profitable (or possible) to move units from s-1,... to p (Phase 3 is done)
                    		xp = nodesPart[p]->UBx - remOfferPartial;//__PHASE4 data
                    		xs = xp;//__PHASE4 data
                    		break;
                    	}
                        currentNodeLL = currentNodeLL->prev;    //will next remove units from a less expensive node (s-1,...)
                    }
					if(currentNodeLL->ptrNode->x == 0){ // this means that the current node has been removed from the LR (not usable anymore with the remaining partials to be tried)
						removeFromLL_backwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->prev
					}
                    else if(currentNodeLL->ptrNode->x < remOfferPartial){ //transfer all the units of currentNodeLL to p. More units still can be transfered from s-1,... to p
                    	LBp2 -= currentNodeLL->ptrNode->x * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
                    	remOfferPartial -= currentNodeLL->ptrNode->x;
                    	if((currentNodeLL->prev) == NULL || (currentNodeLL->prev)->ptrNode->e <= nodesPart[p]->c){ //No longer profitable (or possible) to move units from s-1,... to p (Phase 3 is done)
                    		xs = 0;//__PHASE4 data
                    		xp = nodesPart[p]->UBx - remOfferPartial;//__PHASE4 data
                    		break;
                    	}
                        currentNodeLL = currentNodeLL->prev;    //will next remove units from a less expensive node (s-1,...)
                    }
                    else{ //add remOfferPartial to p to respect UBxp exactly (Phase 3 is done)
                    	LBp2 -= remOfferPartial * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
                    	xs = currentNodeLL->ptrNode->x - remOfferPartial;//__PHASE4 data
                    	remOfferPartial = 0;
                    	xp = nodesPart[p]->UBx;//__PHASE4 data
                    	break;
                    }
                }
            }
        }
		else{                                                //CASE 3) c_p = e_s.
			if(nodesPart[p]->x < nodesPart[p]->LBx){ // A) Must respect the usage bounds. Adjust the LR by transferring units from s,s-1,.. to p (unprofitable or without effect) ->same as 1-A
				//printf("3-A\n");
				int missingUnitsPartial = nodesPart[p]->LBx - nodesPart[p]->x;   //missingUnitsPartial MUST be added to p
				for(;;){
					if(currentNodeLL == NULL){ // Impossible to respect LBxp
						LBp2 = pr->UB+1; // Impossible to respect LBxp
						break;
					}
					if(currentNodeLL->ptrNode->x == 0){ // this means that the current node has been removed from the LR (not usable anymore with the remaining partials to be tried)
						removeFromLL_backwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->prev
					}
                    else if(currentNodeLL->ptrNode->x < missingUnitsPartial){ //transfer all the offer from currentNodeLL to p. Not sufficient yet to respect the bound
                    	LBp2 += currentNodeLL->ptrNode->x * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	missingUnitsPartial -= currentNodeLL->ptrNode->x;
                        currentNodeLL = currentNodeLL->prev;    //will next remove units from a less expensive node (s-1,...)
                    }
                    else{ //add missingUnitsPartial to p to respect LBxp exactly
                    	LBp2 += missingUnitsPartial * (nodesPart[p]->c - currentNodeLL->ptrNode->e);
                    	xs = currentNodeLL->ptrNode->x - missingUnitsPartial;//__PHASE4 data
                    	missingUnitsPartial = 0;
                    	xp = nodesPart[p]->LBx;//__PHASE4 data
                    	break;
                    }
                }
            }
			else if(nodesPart[p]->x > nodesPart[p]->UBx){//B) Must respect the usage bounds. Adjust the LR by transferring units from p to s,s+1,... (unprofitable or without effect) ->same as 2-A
				//printf("3-B\n");
				int exceedingUnitsPartial = nodesPart[p]->x - nodesPart[p]->UBx; //exceedingUnitsPartial MUST be removed from p
				for(;;){
					if(currentNodeLL == NULL){ // Impossible to respect UBxp
						LBp2 = pr->UB+1; 
						break;
					}
					int remOfferCurrent = currentNodeLL->ptrNode->b - currentNodeLL->ptrNode->x; //up to remOfferCurrent units may be transfered from p to currentNodeLL
					if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){// this means that the current node is not usable anymore with the remaining partials to be tried
						removeFromLL_forwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->next
					}
					else if(remOfferCurrent < exceedingUnitsPartial){ //transfer remOfferCurrent units from p to currentNodeLL. More units still MUST be transfered from p to s+1,...
						LBp2 += remOfferCurrent * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
						exceedingUnitsPartial -= remOfferCurrent;
                        currentNodeLL = currentNodeLL->next;    //will next add units on a more expensive node (s+1,...)
                    }
                    else{ //remove exceedingUnitsPartial from p to respect UBxp the bound exactly. (Phase 3 is done)
                    	LBp2 += exceedingUnitsPartial * (currentNodeLL->ptrNode->e - nodesPart[p]->c);
                    	xs = currentNodeLL->ptrNode->x + exceedingUnitsPartial;//__PHASE4 data
                    	exceedingUnitsPartial = 0;//__PHASE4 data
                    	xp = nodesPart[p]->UBx;
                    	break;
                    }
                }
            }
            //printf("3-c\n");
			                                            //C) Usage bounds already respected. No adjustment is required (without effect)
            xp = nodesPart[p]->x;//__PHASE4 data
            xs = currentNodeLL->ptrNode->x;//__PHASE4 data
        }
		//printf("LBp2 = %f\n", LBp2);  //LBp2
		if(LBp2 > nodesPart[p]->LB){  //update LBp 
			nodesPart[p]->LB = LBp2;
		}
		if(nodesPart[p]->LB >= pr->UB){ //p is no longer a partial node candidate
			nodesPart[p]->poss_p = 0;
			pr->nb_poss_p-=1;
			continue; //continue to next node
		}

		// *** Phase 4: adjust the solution from phase 3 to calculate LBxp2 or UBxp2 (depends from the case) ***  O(min(n, b_p/min_{b_j})) for each potential partial node p
		double LBp2Adj = LBp2; // Adjusted LB (on P3|p is partial and ((xp<=xp' (case 1)) or (xp>=xp' (case 2)))
		//special case: p=currentNodeLL, everything removed before p or everyhting added after p in the adjustment phase
		// The adjusted LR from phase 3 is a feasible solution with p partial. 
		if(nodesPart[p]->i == currentNodeLL->ptrNode->i || currentNodeLL == NULL){
			nodesPart[p]->LBx = xp; //this is exactly the optimal solution to (P3|p partial)
			nodesPart[p]->UBx = xp;
    		pr->UB = LBp2; //Adjustment of UB (incumbent solution)
    		//TODO: pr->incumbent = z; 
    		printf("%New incumbent solution: feasible solution to P3 from LR with partial node %d\n", nodesPart[p]->i);
    	}

		//general case 1: p!=currentNodeLL and c_p < e_currentNodeLL (we move units from p to currentNodeLL, currentNodeLL+1... until LBp' reaches UB or xp' reaches LBxp) -> LBxp2
    	else if(nodesPart[p]->c < currentNodeLL->ptrNode->e){
    		int LBxp2;
    		//move units from p to currentNodeLL (s)
    		double cost_by_unit_moved = currentNodeLL->ptrNode->e - nodesPart[p]->c;
    		int max_units = floor((pr->UB-LBp2Adj)/cost_by_unit_moved);//maximum number of units that can be moved from p to s without exceeding the UB
    		int remOfferSplit = currentNodeLL->ptrNode->b-xs;
    		if(max_units <= remOfferSplit){ //moving more than max_units from p to s leads to exceed the current UB and s has >= max_units remaining offer left on it. We have the UBxp 
    			LBxp2 = xp - max_units;
    			if(LBxp2 > nodesPart[p]->LBx){ //improvement of LBxp
    				nodesPart[p]->LBx = LBxp2;
    			}
    			continue; //continue with next potential partial node 
    		}
    		else{ // we move remOfferSplit units from p to s (fill s) and continue on the next node (s+1)
    			LBp2Adj += cost_by_unit_moved*remOfferSplit;
    			xp -= remOfferSplit;
    			currentNodeLL = currentNodeLL->next;
    		}
    		for(;;){
    			if(xp <= nodesPart[p]->LBx){ //We did not improved LBxp
    				break;
    			}
    			if(currentNodeLL == NULL){ //No longer possible to move units from p to s+1,... (Phase 4 is done)
    				LBxp2 = xp;
    				if(LBxp2 > nodesPart[p]->LBx){ //improvement of LBxp
    					nodesPart[p]->LBx = LBxp2;
    				}
    				break;
    			}
    			if(nodesPart[p]->i == currentNodeLL->ptrNode->i){ //skip currentNodeLL if it is the same node as p
    				currentNodeLL = currentNodeLL->next;
    			}
    			else if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){ // this means that the current node has been removed from the LR (not usable anymore with the remaining partials to be tried)
					removeFromLL_forwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->next
				}
				else{
					cost_by_unit_moved = currentNodeLL->ptrNode->e - nodesPart[p]->c;
		    		max_units = floor((pr->UB-LBp2Adj)/cost_by_unit_moved);//maximum number of units that can be moved from p to currentNodeLL without exceeding the UB
		    		if(max_units <= currentNodeLL->ptrNode->b){ //moving more than max_units from p to currentNodeLL leads to exceed the current UB and currentNodeLL has >= max_units offer. We have the UBxp 
		    			LBxp2 = xp - max_units;
    					if(LBxp2 > nodesPart[p]->LBx){ //improvement of LBxp
    						nodesPart[p]->LBx = LBxp2;
    					}
    					break;
    				}
    				else{ // we move currentNodeLL->ptrNode->b units from p to currentNodeLL (fill currentNodeLL) and continue on the next node (currentNodeLL+1)
    					LBp2Adj += cost_by_unit_moved*currentNodeLL->ptrNode->b;
    					xp -= currentNodeLL->ptrNode->b;
    					currentNodeLL = currentNodeLL->next;
    				}
    			}

    		}
    	}
		//general case 2: p!=currentNodeLL and c_p > e_currentNodeLL (we move units from currentNodeLL, currentNodeLL-1... to p until LBp' reaches UB or xp' reaches UBxp) -> UBxp2
    	else{
    		int UBxp2;
    		//move units from currentNodeLL (s) to p
    		double cost_by_unit_moved = nodesPart[p]->c - currentNodeLL->ptrNode->e;
    		int max_units = floor((pr->UB-LBp2Adj)/cost_by_unit_moved);//maximum number of units that can be moved from s to p without exceeding the UB
    		if(cost_by_unit_moved == 0){ // Division by 0 above. cost_by_unit_moved being 0, as many units as we want can be moved without exceeding the UB
    			max_units = pr->D;
    		}
    		if(max_units <= xs){ //moving more than max_units from s to p leads to exceed the current UB and s has >= max_units used on it. We have the UBxp 
    			UBxp2 = xp + max_units;
    			if(UBxp2 < nodesPart[p]->UBx){ //improvement of UBxp
    				nodesPart[p]->UBx = UBxp2;
    			}
    			continue; //continue with next potential partial node 
    		}
    		else{ // we move xs units from s to p and continue on the previous node (s-1)
    			LBp2Adj += cost_by_unit_moved*xs;
    			xp += xs;
    			currentNodeLL = currentNodeLL->prev;
    		}
    		for(;;){
    			if(xp >= nodesPart[p]->UBx){ //We did not improved UBxp
    				break;
    			}
				if(currentNodeLL == NULL){ //No longer possible to move units from s-1,... to p (Phase 4 is done)
					UBxp2 = xp;
    				if(UBxp2 < nodesPart[p]->UBx){ //improvement of UBxp
    					nodesPart[p]->UBx = UBxp2;
    				}
    				break;
    			}
    			cost_by_unit_moved = nodesPart[p]->c - currentNodeLL->ptrNode->e;
				max_units = floor((pr->UB-LBp2Adj)/cost_by_unit_moved);//maximum number of units that can be moved from s to p without exceeding the UB
				if(cost_by_unit_moved == 0){ // Division by 0 above. cost_by_unit_moved being 0, as many units as we want can be moved without exceeding the UB
	    			max_units = pr->D;
	    		}
				if(partial_node_comparator(nodesPart[p], currentNodeLL->ptrNode) < 0){ // this means that the current node has been removed from the LR (not usable anymore with the remaining partials to be tried)
					removeFromLL_backwards(&currentNodeLL);//currentNodeLL is removed from the linked list and currentNodeLL becomes currentNodeLL->prev
				}
				else if(max_units <= currentNodeLL->ptrNode->b){//moving more than max_units from currentNodeLL to p leads to exceed the current UB and currentNodeLL has >= max_units used on it. We have the UBxp 
					UBxp2 = xp + max_units;
    				if(UBxp2 < nodesPart[p]->UBx){ //improvement of UBxp
    					nodesPart[p]->UBx = UBxp2;
    				}
    				break;
    			}
    			else{ // we move all the units from currentNodeLL to p and continue on the previous node (currentNodeLL-1) 
    				LBp2Adj += cost_by_unit_moved*currentNodeLL->ptrNode->b;
    				xp += currentNodeLL->ptrNode->b;
    				currentNodeLL = currentNodeLL->prev;
    			}
    		}
    	}
    }

	//free(*nodesLR);   // free the heap
	//free(*nodesPart); // free the arraylist 
	util_free_LL(firstNodeLL, lastNodeLL); // free the linkedlist
}

//=================================== 10.3-10.4: Bounds on (xp|p is partial) and dominance relations ===================================
    void LB_UBxp1(int p, Node * nodes, P3struct * P3L, Prob * pr, int sizeP3L){
    	double bound;
    	for (int i = 0; i < sizeP3L; ++i)
    	{
    		if(P3L[i].surplus>pr->ub) continue;
		if(P3L[i].P3lambda == 0) return; //no more value in L
		if(P3L[i].lambda>nodes[p].c){
			bound = (int)nodes[p].b - floor((P3L[i].P3lambda - pr->UB)/(nodes[p].c - P3L[i].lambda));
			if(bound>nodes[p].LBx){
				nodes[p].LBx = bound;
				if(nodes[p].LBx > nodes[p].UBx){
					nodes[p].poss_p = 0;
					pr->nb_poss_p-=1;
					return;
				}
			} 
		}
		if(P3L[i].lambda<nodes[p].c){
			//printf("%f\n", P3L[i].P3lambda - (double)pr->UB);
			bound = nodes[p].b - ceil((P3L[i].P3lambda - pr->UB)/(nodes[p].c - P3L[i].lambda));
			if(bound<nodes[p].UBx){
				nodes[p].UBx = bound;
				if(nodes[p].LBx > nodes[p].UBx){
					nodes[p].poss_p = 0;
					pr->nb_poss_p-=1;
					return;
				}
			} 
		}
	}
}

//Dominance relation. Faster implementation using LinkedList. Only use potential partial nodes as potential dominant partial nodes
void LB_UBxp1_all(Node * nodes, P3struct * P3L, Prob * pr, int n, int sizeP3L){
	if(pr->nb_poss_p == 0)return;
	for (int i = 0; i < n; ++i)
	{
		if(nodes[i].poss_p) LB_UBxp1(i, nodes, P3L, pr, sizeP3L);
	}
}

//Dominance relation. Faster implementation using LinkedList. Only use potential partial nodes as potential dominant partial nodes
void LBxp3_all_efficient(Node * nodes, int n, Prob * pr){ 
	if(pr->nb_poss_p == 0)return;

    // nodesPart: Array list of pointers to nodes that may be partial in an optimal solution to P3, sorted in 
    //			  increasing order (using total order from the article)
	Node *nodesPart[n];
	int nb_poss_part;
    build_arraylist_partials(n, nodes, nodesPart, &nb_poss_part); //O(n logn)

    DominanceNodeLinkedList *firstDomNodeLL, *currentDomNodeLL;
	firstDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // first dominanceNodeLL

	firstDomNodeLL->ptrNode = nodesPart[nb_poss_part-1]; // the potential partial node with the highest c_j first dominates everywhere it is valid (for any use x_p between 0 and b_j)
	firstDomNodeLL->xMax = nodesPart[nb_poss_part-1]->b; // up to b_p units can be used on this node
	firstDomNodeLL->next = NULL;

	double xCritical;
	for (int p = nb_poss_part-2; p >= 0; --p){
		if(nodesPart[p]->UBx >= nodesPart[p]->b){
			nodesPart[p]->UBx = nodesPart[p]->b-1; //by definition, a partial node cannot be completely used. Hence, nodesPart[p]->UBx  is set to b_p-1
		}

    	//printf(" p = %d\n", nodesPart[p]->i);
    	int LBxp3 = 1;                                     // reset LBxp3 for the new node p
    	if(nodesPart[p]->f <= firstDomNodeLL->ptrNode->f){ // p dominates for any use in [0,b_p]. We must update firstDomNodeLL

    		DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
    		newDomNodeLL->ptrNode = nodesPart[p];
    		newDomNodeLL->xMax = nodesPart[p]->b;

    		DominanceNodeLinkedList *nextDomNodeLL = firstDomNodeLL; //free memory on the destroyed LLnodes (all the previously dominant nodes with b_j <= b_p)
    		while(nextDomNodeLL != NULL && nextDomNodeLL->xMax<=nodesPart[p]->b){ //nextDomNode can be destroyed (is dominated by p on all its dominance interval)
    			DominanceNodeLinkedList *temp = nextDomNodeLL->next;
    			free(nextDomNodeLL);        
    			nextDomNodeLL = temp;
    		}

			newDomNodeLL->next = nextDomNodeLL; // nextDomNodeLL=NULL if b_p >= b_j for all j (previously dominant nodes on a given sub-interval)
			firstDomNodeLL = newDomNodeLL;
		}
		else{
			currentDomNodeLL = firstDomNodeLL;
    		//while(nodesPart[p]->c < currentDomNodeLL->ptrNode->c){ //2nd condition is always true in general. This check is just to avoid a division by 0 if two nodes have exactly the same c_j
			for(;;){

				xCritical = (nodesPart[p]->f-currentDomNodeLL->ptrNode->f)/(currentDomNodeLL->ptrNode->c-nodesPart[p]->c);

    			if(xCritical <= nodesPart[p]->b && xCritical <= currentDomNodeLL->xMax){ //Case 1) p is dominated on [0,xCritical] and dominates on [xCritical,b_p]
    				// current[..,xCritical] -> new[xCritical,b_p] -> next[b_p,..]
    				DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
    				DominanceNodeLinkedList * nextDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // next dominanceNodeLL (dominates after b_p)
    				
    				newDomNodeLL->ptrNode = nodesPart[p];
    				newDomNodeLL->xMax = nodesPart[p]->b;
    				
    				nextDomNodeLL->ptrNode = currentDomNodeLL->ptrNode;
    				nextDomNodeLL->xMax = currentDomNodeLL->xMax;
    				nextDomNodeLL->next = currentDomNodeLL->next;
    				while(nextDomNodeLL!=NULL && nextDomNodeLL->xMax <= nodesPart[p]->b){
    					DominanceNodeLinkedList *temp = nextDomNodeLL->next;
    					free(nextDomNodeLL);        
    					nextDomNodeLL = temp;
    				}
    				
    				newDomNodeLL->next = nextDomNodeLL;
    				currentDomNodeLL->next = newDomNodeLL;
    				currentDomNodeLL->xMax = xCritical;

    				LBxp3 = ceil(xCritical);
    				break;
    			}

    			else if(xCritical <= nodesPart[p]->b && xCritical >= currentDomNodeLL->ptrNode->b){ //current dominates p on [0,b_current]. Must continue for [b_current,b_p]
    				LBxp3 = currentDomNodeLL->ptrNode->b;

    				if(currentDomNodeLL->next == NULL){
                        DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
                        DominanceNodeLinkedList * nextDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // next dominanceNodeLL (dominates after b_p)

                        newDomNodeLL->ptrNode = nodesPart[p];
                        newDomNodeLL->xMax = nodesPart[p]->b;
                        newDomNodeLL->next = NULL;
                        currentDomNodeLL->next = newDomNodeLL;
                        break;
                    }
                    else if(currentDomNodeLL->next->ptrNode->f >= nodesPart[p]->f){
    					DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
    					DominanceNodeLinkedList * nextDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // next dominanceNodeLL (dominates after b_p)

    					newDomNodeLL->ptrNode = nodesPart[p];
    					newDomNodeLL->xMax = nodesPart[p]->b;

    					nextDomNodeLL->ptrNode = currentDomNodeLL->next->ptrNode;
    					nextDomNodeLL->xMax = currentDomNodeLL->next->xMax;
    					nextDomNodeLL->next = currentDomNodeLL->next->next;
    					while(nextDomNodeLL!=NULL && nextDomNodeLL->xMax <= nodesPart[p]->b){
    						DominanceNodeLinkedList *temp = nextDomNodeLL->next;
    						free(nextDomNodeLL);        
    						nextDomNodeLL = temp;
    					}

    					newDomNodeLL->next = nextDomNodeLL;
    					currentDomNodeLL->next = newDomNodeLL;
    					break;
    				}
    				currentDomNodeLL=currentDomNodeLL->next;
    				//printf("a\n");
    			}

    			else if(xCritical <= nodesPart[p]->b){ // xCritical <= b_p &&  xMax < xCritical < b_current
    				LBxp3 = ceil(xCritical);
    				currentDomNodeLL=currentDomNodeLL->next;
    				//printf("b\n");
    			}

    			else if(xCritical >= nodesPart[p]->b && xCritical <= currentDomNodeLL->ptrNode->b){ //current dominates p on [0,b_p]. End
    				LBxp3 = ceil(xCritical);
    				break;
    			}

    			else{ //(xCritical >= nodesPart[p]->b && xCritical >= currentDomNodeLL->ptrNode->b //current dominates p on [0,b_current].
    				if(currentDomNodeLL->ptrNode->b >= nodesPart[p]->b){ //current dominates p everywhere
    					LBxp3 = currentDomNodeLL->ptrNode->b;
    					break;
    				}
    				else{ //current dominates p on [0,b_current]. Must continue for [b_current,b_p]
    					LBxp3 = currentDomNodeLL->ptrNode->b;
    					if(currentDomNodeLL->next == NULL){
	                        DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
	                        DominanceNodeLinkedList * nextDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // next dominanceNodeLL (dominates after b_p)

	                        newDomNodeLL->ptrNode = nodesPart[p];
	                        newDomNodeLL->xMax = nodesPart[p]->b;
	                        newDomNodeLL->next = NULL;
	                        currentDomNodeLL->next = newDomNodeLL;
	                        break;
	                    }

	                    else if(currentDomNodeLL->next->ptrNode->f >= nodesPart[p]->f){
    					DominanceNodeLinkedList * newDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // new dominanceNodeLL (p dominates after xCritical)
    					DominanceNodeLinkedList * nextDomNodeLL = (DominanceNodeLinkedList*)malloc(sizeof(DominanceNodeLinkedList)); // next dominanceNodeLL (dominates after b_p)

    					newDomNodeLL->ptrNode = nodesPart[p];
    					newDomNodeLL->xMax = nodesPart[p]->b;

    					nextDomNodeLL->ptrNode = currentDomNodeLL->next->ptrNode;
    					nextDomNodeLL->xMax = currentDomNodeLL->next->xMax;
    					nextDomNodeLL->next = currentDomNodeLL->next->next;
    					while(nextDomNodeLL!=NULL && nextDomNodeLL->xMax <= nodesPart[p]->b){
    						DominanceNodeLinkedList *temp = nextDomNodeLL->next;
    						free(nextDomNodeLL);        
    						nextDomNodeLL = temp;
    					}

    					newDomNodeLL->next = nextDomNodeLL;
    					currentDomNodeLL->next = newDomNodeLL;
    					break;
    				}
    				currentDomNodeLL=currentDomNodeLL->next;
    					//printf("c\n");
    			}
    		}
    	}
    }
    if(LBxp3 >= nodesPart[p]->LBx){
    	nodesPart[p]->LBx = LBxp3;
    	if(nodesPart[p]->LBx > nodesPart[p]->UBx){
    		nodesPart[p]->poss_p = 0;
    		pr->nb_poss_p-=1;
    	}
    }
}

while(firstDomNodeLL != NULL){ //clean memory
	DominanceNodeLinkedList *temp = firstDomNodeLL->next;
	free(firstDomNodeLL);        
	firstDomNodeLL = temp;
}
}

//Apply dominance relation on a node p (inefficient implementation)
static void LBxp3(int p, Node * nodes, int n, Prob * pr){ //O(n) for each node
	int64_t lb;
	for (int q = 0; q < n; ++q)
	{
		if(nodes[q].poss_p == 0 || nodes[q].c <= nodes[p].c || nodes[q].f >= nodes[p].f) continue; //q can't dominate p in these cases
		lb = ceil((nodes[p].f - nodes[q].f)/(nodes[q].c - nodes[p].c));
		if(nodes[q].b < lb){
			lb = nodes[q].b;
		}
		if(lb > nodes[p].LBx){
			nodes[p].LBx=lb;
			if(nodes[p].LBx > nodes[p].UBx){
				nodes[p].poss_p = 0;
				pr->nb_poss_p-=1;
				break;
			}
		}
	}
}

//Dominance relation for partial nodes. O(n^2) implementation. Only tests the dominance between potential partial nodes
void LBxp3_all(Node * nodes, int n, Prob * pr){ 
	if(pr->nb_poss_p == 0)return;

	for (int i = 0; i < n; ++i)
		if(nodes[i].poss_p) LBxp3(i, nodes, n, pr);
}

//======================================================= Simple functions =======================================================

// Reject p as a partial node if cp is not in (Cmin, Cmax). 
// Update LB = min{LBp | p is potentially partial}
void CminCmaxLBFilter_LBUpdate(Prob *pr, Node *nodes){
	double min_LBp = pr->UB;
	double min_cp  = pr->uc;
	double max_cp  = pr->lc;
	for (int i = 0; i < pr->n; ++i)
	{
		if(!nodes[i].poss_p) continue;
		if(nodes[i].c<pr->Cmin || nodes[i].c>pr->Cmax || (nodes[i].c==pr->Cmin && pr->cmin_closed_interval==0) || (nodes[i].c==pr->Cmax && pr->cmax_closed_interval==0)){
			nodes[i].poss_p=0;
			pr->nb_poss_p-=1;
			continue;
		}

		min_LBp = MIN(min_LBp, nodes[i].LB);
		min_cp = MIN(min_cp, nodes[i].c);
		max_cp = MAX(max_cp, nodes[i].c);
	}

	pr->LB = MAX(pr->LB, min_LBp);

	if(min_cp > pr->Cmin){
		pr->Cmin = min_cp;
		pr->cmin_closed_interval = 1; //this node can stricly improve the incumbent solution as a partial
	}

	if(max_cp < pr->Cmax){
		pr->Cmax = max_cp;
		pr->cmax_closed_interval = 1; //this node can stricly improve the incumbent solution as a partial
	}
}

// Computes a slightly modified version of Dantzig's upper bound (sum_{i=1}^s u_j - (c_p-c_s)*(min(b_s-x_s, b_p))),
// where s is the split node and p := argmin_{i<=s}(c_i).
void LB_UB_Dantzig_modif(Prob *pr, Node *nodes, int computeSol){
	int demand = pr->D;
	double dantzigUB = 0;
	int p = 0;
	for (int i = 0; i < pr->n; ++i)
	{
		if(demand > nodes[i].b){
			dantzigUB+=nodes[i].u;
			demand -= nodes[i].b;
			if(nodes[i].c > nodes[p].c){
				p = i;
			}
		}
		else{
			pr->s = i;
			pr->xs = demand;
			pr->LB_linear = dantzigUB + demand*nodes[i].e;
			dantzigUB += nodes[i].f+demand*nodes[i].c;
			if(nodes[i].c <= nodes[p].c){
				dantzigUB -= MIN(nodes[i].b-demand,nodes[p].b) * (nodes[p].c - nodes[i].c);
			}
			break;
		}
	}

	pr->LB = MAX(pr->LB, pr->LB_linear);
	if(dantzigUB < pr->UB){
		pr->UB = dantzigUB;
		//=== Update the incumbent solution ===
		if(computeSol){
			int offer_comp = 0;
			for (int i = 0; i <= pr->s; ++i){ 
				if(i != p){
					int x_i = nodes[i].b;
					offer_comp += x_i;
					pr->incumbent[i] = x_i; 
				}
			}
			pr->incumbent[p] = pr->D - offer_comp;
		}
		//======================================
	}
	//printf("UB/LB = %f\n", pr->UB/pr->LB);
	//printf("epgap=10^{-6}=0.000001 LB+UB*epgap > UB ? : %d \n",(pr->LB+0.000001*pr->UB > pr->UB));
}

// Simple linear filtering in O(n) (given that the nodes are already sorted in non-decreasing order of linearized cost e_j)
// A node j < s can't be unused          (poss_0 = FALSE) if linearLB+b_j(e_s - e_j) >= UB
// A node j > s can't be completely used (poss_1 = FALSE) if linearLB+b_j(e_j - e_s) >= UB
void linear_filtering(Prob *pr, Node *nodes){
	int s = pr->s;
	int nb_poss_0_FALSE = 0;
	int nb_poss_1_FALSE = 0;
	for (int i = 0; i < s; ++i)
	{
		if(nodes[i].poss_0 == 0){ //we already know that i can't be unused
			continue;
		}
		if(pr->LB_linear + nodes[i].b*(nodes[s].e-nodes[i].e) > pr->UB){
			nodes[i].poss_0 = 0;
			pr->nb_poss_0 -=1;
		}
	}
	for (int i = s+1; i < pr->n; ++i)
	{
		if(nodes[i].poss_1 == 0){ //we already know that i can't be completely used
			continue;
		}
		if(pr->LB_linear + nodes[i].b*(nodes[i].e-nodes[s].e) > pr->UB){
			nodes[i].poss_1 = 0;
			pr->nb_poss_1 -=1;
		}
	}
}

// Linear filtering in O(n * min(n,max{b_j}/min{b_j})) (given that the nodes are already sorted in non-decreasing order of linearized cost e_j)
// A node j < s can't be unused          (poss_0 = FALSE) if linearLB - u_j + sum_{k=s}^t(e_k * m_k) >= UB, where m_k is the number of units that are moved from j to k. x_s+m_s<=b_s and m_k<=b_k for k \in {s+1,s+2,...,t}
// A node j > s can't be completely used (poss_1 = FALSE) if linearLB + u_j - sum_{k=s:-1}^t(e_k * r_k) >= UB, where r_k is the number of units that are removed from k. r_s<=x_s and r_k<=b_k for k \in {s-1,,s-2,...,t}
void linear_filtering_improved(Prob *pr, Node *nodes){
	int s = pr->s;
	int xs = pr->xs;
	int n = pr->n;
	int nb_poss_0_FALSE = 0;
	int nb_poss_1_FALSE = 0;
	int remaining_offer_s = nodes[s].b-xs;
	for (int i = 0; i < s; ++i) //fixing nodes[i].poss_0 to FALSE for nodes j<s that can't be unused without the linear LB to excess the UB
	{
		if(nodes[i].poss_0 == 0){ //we already know that i can't be unused
			continue;
		}
		double LB_i_unused = pr->LB_linear-nodes[i].u;
		int xi = nodes[i].b;
		if(xi <= remaining_offer_s){
			LB_i_unused += xi*nodes[s].e;
			xi = 0;
		}
		else{
			LB_i_unused += remaining_offer_s*nodes[s].e;
			xi -= remaining_offer_s;
			int k = s+1;
			while(k<n && xi > nodes[k].b){
				LB_i_unused += nodes[k].b*nodes[k].e;
				xi -=  nodes[k].b;
				k++;
			}
			if(k<n){
				LB_i_unused += xi*nodes[k].e;
				xi = 0;
			}
		}
		if(xi>0 || LB_i_unused > pr->UB){
			nodes[i].poss_0 = 0;
			pr->nb_poss_0 -=1;
		}
	}

	for (int i = s+1; i < n; ++i) //fixing nodes[i].poss_1 to FALSE for nodes j>s that can't be completely used without the linear LB to excess the UB
	{
		if(nodes[i].poss_1 == 0){ //we already know that i can't be completely used
			continue;
		}
		double LB_i_completely_used = pr->LB_linear+nodes[i].u;
		int exceeding_offer = nodes[i].b;
		if(exceeding_offer <= xs){
			LB_i_completely_used -= exceeding_offer*nodes[s].e;
			exceeding_offer = 0;
		}
		else{
			LB_i_completely_used -= xs*nodes[s].e;
			exceeding_offer -= xs;
			int k = s-1;
			while(k > 0 && exceeding_offer > nodes[k].b){
				LB_i_completely_used -= nodes[k].b*nodes[k].e;
				exceeding_offer -=  nodes[k].b;
				k--;
			}
			if(k > 0){
				LB_i_completely_used -= exceeding_offer*nodes[k].e;
				exceeding_offer = 0;
			}
		}
		if(exceeding_offer > 0 || LB_i_completely_used > pr->UB){
			nodes[i].poss_1 = 0;
			pr->nb_poss_1 -=1;
		}
	}
}

// Linear filtering in O(n * min(n,max{b_j}/min{b_j})) (given that the nodes are already sorted in non-decreasing order of linearized cost e_j)
// A node j < s can't be unused          (poss_0 = FALSE) if linearLB - u_j + sum_{k=s}^t(e_k * m_k) >= UB, where m_k is the number of units that are moved from j to k. x_s+m_s<=b_s and m_k<=b_k for k \in {s+1,s+2,...,t}
// A node j > s can't be completely used (poss_1 = FALSE) if linearLB + u_j - sum_{k=s:-1}^t(e_k * r_k) >= UB, where r_k is the number of units that are removed from k. r_s<=x_s and r_k<=b_k for k \in {s-1,,s-2,...,t}
// ONLY CONSIDER NODES j \prec q to improve filtering linear_filtering_improved a bit more
void linear_filtering_improved_Q(Prob *pr, Node *nodes){ 
	int n = pr->n;
	int q = -1;
	int sq = 0;
	int xsq = 0;
	double Zlp_q = 0;
	int inQ[n];  //inQ[i]==1 iff node i is in set Q

	// find the maximal element q among partial nodes candidates
	for (int i = 0; i < n; ++i)
	{
		if(nodes[i].poss_p){
			if ( q==-1 || i_prec_j(q, i, nodes) ){
				q = i; 
			}
		}
	}

	// calculate the linear relaxation on the subset of nodes Q corresponding to the subset of nodes j such that j<=q according to the total order
	int demand = pr->D;
	for (int i = 0; i < n; ++i)
	{
		if(i_prec_j(q, i, nodes)){ // q < i, hence i cannot be used in a solution with a partial node p \in P 
			inQ[i] = 0;
			nodes[i].poss_1 == 0; // i cannot be complete
			continue;
		}
		inQ[i] = 1;

		if(demand > nodes[i].b){
			Zlp_q += nodes[i].u;
			demand -= nodes[i].b;
		}
		else{
			for (int j = sq+1; j < n; ++j) // finish to construct set Q by checking nodes j \in {sq+1,...,n-1}
			{
				if(i_prec_j(q, j, nodes)){
					inQ[j] = 0;
					nodes[j].poss_1 == 0; // j cannot be complete
				}
				else{
					inQ[j] = 1;
				}
			}
			sq = i;
			xsq = demand;
			Zlp_q += demand *nodes[sq].e;
			demand -= demand;
			break;
		}
	}

	if(Zlp_q >= pr->UB || demand > 0){
		printf("*****************************************************************************************************\n");
		printf("* linear_filtering_improved2: NO SOLUTION USING ONLY NODES ON SET Q CAN IMPROVE THE OPTIMAL SOLUTION*\n");
		printf("*****************************************************************************************************\n");
		printf("Zlp_q >= pr->UB : %d\n",Zlp_q >= pr->UB );
		printf("demand > 0 : %d\n",demand > 0 );
		return;
	}


	int s = sq;
	int xs = xsq;
	int remaining_offer_s = nodes[s].b-xs;
	for (int i = 0; i < s; ++i) //fixing nodes[i].poss_0 to FALSE for nodes j<s that can't be unused without the linear LB to excess the UB
	{
		if(nodes[i].poss_0 == 0){ //we already know that i can't be unused
			continue;
		}
		double LB_i_unused = pr->LB_linear-nodes[i].u;
		int xi = nodes[i].b;
		if(xi <= remaining_offer_s){
			LB_i_unused += xi*nodes[s].e;
			xi = 0;
		}
		else{
			LB_i_unused += remaining_offer_s*nodes[s].e;
			xi -= remaining_offer_s;
			int k = s+1;
			while(k<n && xi > nodes[k].b){
				if(inQ[i] == 1){  //otherwise, j in not in Q. We ignore it
					LB_i_unused += nodes[k].b*nodes[k].e;
					xi -=  nodes[k].b;
				}
				k++;
			}
			if(k<n){
				LB_i_unused += xi*nodes[k].e;
				xi = 0;
			}
		}
		if(xi>0 || LB_i_unused > pr->UB){ // no feasible solution using only nodes on set Q in which i is unused, or no solution that could improve the incumbent
			nodes[i].poss_0 = 0;
			pr->nb_poss_0 -=1;
		}
	}

	for (int i = s+1; i < n; ++i) //fixing nodes[i].poss_1 to FALSE for nodes j>s that can't be completely used without the linear LB to excess the UB
	{
		if(nodes[i].poss_1 == 0){ //we already know that i can't be completely used
			continue;
		}
		double LB_i_completely_used = pr->LB_linear+nodes[i].u;
		int exceeding_offer = nodes[i].b;
		if(exceeding_offer <= xs){
			LB_i_completely_used -= exceeding_offer*nodes[s].e;
			exceeding_offer = 0;
		}
		else{
			LB_i_completely_used -= xs*nodes[s].e;
			exceeding_offer -= xs;
			int k = s-1;
			while(k > 0 && exceeding_offer > nodes[k].b){
				if(inQ[i] == 1){   //otherwise, j in not in Q. It has not been used in the linear relaxation, hence we cannot remove units from it
					LB_i_completely_used -= nodes[k].b*nodes[k].e;
					exceeding_offer -=  nodes[k].b;
				}
				k--;
			}
			if(k > 0){
				LB_i_completely_used -= exceeding_offer*nodes[k].e;
				exceeding_offer = 0;
			}
		}
		if(exceeding_offer > 0 || LB_i_completely_used > pr->UB){ // no feasible solution using only nodes on set Q in which i is completely used, or no solution that could improve the incumbent
			nodes[i].poss_1 = 0;
			pr->nb_poss_1 -=1;
		}
	}

}

//======================================================= Output functions =======================================================
//print basic information about node i
void printNode(Node * nodes, int i){
	printf("%s%d\n", "i = ",nodes[i].i);
	printf("%s%d\n", "b = ", nodes[i].b);
	printf("%s%f\n", "c = ",nodes[i].c);
	printf("%s%f\n", "f = ", (double)nodes[i].f);
	printf("%s%f\n","e = ", nodes[i].e);
	printf("%s%f\n","u = ", nodes[i].u);
	printf("%s%f\n","LB = ", nodes[i].LB);
	printf("%s%d\n", "LBx = ",nodes[i].LBx);
	printf("%s%d\n\n", "UBx = ",nodes[i].UBx);
}

//print detailed information about node i
void printNodeDetailed(Node * nodes, int i){
	printf("%s%d\n", "i = ",nodes[i].i);
	printf("%s%d\n", "b = ", nodes[i].b);
	printf("%s%f\n", "c = ",nodes[i].c);
	printf("%s%f\n", "f = ", nodes[i].f);
	printf("%s%f\n","e = ", nodes[i].e);
	printf("%s%f\n","u = ", nodes[i].u);
	printf("%s%f\n","LB = ", nodes[i].LB);
	printf("%s%d\n", "LBx = ",nodes[i].LBx);
	printf("%s%d\n", "UBx = ",nodes[i].UBx);
	printf("%s%d\n", "x = ", nodes[i].x);
	printf("%s%d\n\n", "poss_p = ", nodes[i].poss_p);
}

//print basic information about all the nodes
void printNodes(Node * nodes, int n){
	for (int i = 0; i < n; ++i)
	{
		printNode(nodes, i);
	}
}

//print detailed information about all the nodes
void printNodesDetailed(Node * nodes, int n){
	for (int i = 0; i < n; ++i)
	{
		printNodeDetailed(nodes, i);
	}
}

//print important information about the instance
void printProblem(Prob *pr){
	char cmin_bound = (pr->cmin_closed_interval)?'[':'(';
	char cmax_bound = (pr->cmax_closed_interval)?']':')';
	printf("%s\n",   "============================" );
	printf("%s %d\n","n         = ",pr->n);
	printf("%s %d\n","D         = ",pr->D);
	printf("%s %c%f%s\n","Cmin      =",cmin_bound ,pr->Cmin,",");  
	printf("%s %f%c\n","Cmax      = ",pr->Cmax,cmax_bound); 
	printf("%s %f\n","LB        = ",pr->LB);
	if(pr->UB <  DBL_MAX){
		printf("%s %f\n","UB        = ",pr->UB);
	}
	else{
		printf("%s\n","UB        =  Infinity");
	}
	printf("%s %d\n","nb_poss_p = ",pr->nb_poss_p);
	printf("%s %d\n","nb_poss_1 = ",pr->nb_poss_1);
	printf("%s %d\n","nb_poss_0 = ",pr->nb_poss_0);
	printf("%s\n",   "============================" );
}

//print the incumbent solution
void printSolution(Prob *pr, Node *nodes){
	printf("Optimal solution\n----------------------------");
	for (int i = 0; i < pr->n; ++i)
	{
		printf("x[%d]* = %d \n", i, pr->incumbent[i]);
	}
}

//compute and return the excess supply on the nodes that are used in the optimal solution
int getExcessSupply(Prob *pr, Node *nodes){
	int offer_comp = 0; // total offer on used nodes (including the partial)
	for (int i = 0; i < pr->n; ++i)
	{
		if(pr->incumbent[i] > 0){
			offer_comp += nodes[i].b;
		}
	}
	return (offer_comp - pr->D);
}

//print the nodes that may improve the incumbent solution while being partials
void printPotentialPartials(Node * nodes, int n){
	for (int i = 0; i < n; ++i)
	{
		if(nodes[i].poss_p)
			printNode(nodes, i);
	}
}

//output the information related to struct P3L
void printP3L(P3struct * P3L, int sizeP3L){
	printf("%s\n",   "----------------------------" );
	for (int i = 0; i < sizeP3L; ++i)
	{
		if(P3L[i].lambda == 0) break;
		printf("%s %f\n","lambda      = ",P3L[i].lambda);
		printf("%s %f\n","cp          = ",P3L[i].cp);
		printf("%s %d\n","p           = ",P3L[i].p);
		printf("%s %f\n","P3(lambda)  = ",P3L[i].P3lambda);
		printf("%s %f\n","Zssfctp(z*) = ",P3L[i].Zssfctp);
		printf("%s %d\n","surplus     = ",P3L[i].surplus);
		printf("%s\n",   "----------------------------" );
	}
}

//output the information related to struct LR
void printLR(LRstruct * LR){
	printf("%s\n",   "----------------------------" );
	printf("%s %f\n","cost   = ",LR->cost);
	printf("%s %d\n","s      = ",LR->s);
	printf("%s %d\n","xs     = ",LR->xs);
	printf("%s\n",   "----------------------------" );
}


//output to be used to generate the same problem in Julia
void outputJulia(Node * nodes, int n){
	printf("%s\n", "Nodes=Vector{Node}()");
	for (int i = 0; i < n; i++ ) 
	{
		printf("push!(Nodes,Node(%d",i+1);
		printf(",%d",nodes[i].b);
		printf(",%f",nodes[i].c);
		printf(",%f",nodes[i].f);
		printf(",%f",nodes[i].e);
		printf(",%f",nodes[i].u);
		printf(",0,0,0,0,0,0))\n");
	}
} 

//write the problem in a file that is usable with Klose's code
void outputFileKlose(Node * nodes, int n, int D, int file_idx){
	/* File pointer to hold reference to our file */
	FILE * outfile;
	char pname[1024]= "./output_Klose/"; 
	char str[1024];
	strcat( pname, "outfile" );
	sprintf(str, "%d", file_idx);
	strcat( pname, str ); 
	strcat( pname, ".fctp" );
	//printf("Generating file %s\n",pname);
	outfile= fopen( pname, "wt" );
	fprintf(outfile, "%d\n", n );
	fprintf(outfile, "%d\n", D );
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%d ",nodes[i].b );
	}
	fprintf(outfile,"\n");

	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%15.16f ", nodes[i].c );
	}
	fprintf(outfile,"\n");

	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%d ",nodes[i].f );
	}
	fclose( outfile );
}

//write the problem in an .mps file that is usable by Gurobi
void outputFileGurobi(Node * nodes, int n, int D, char * gen_str, int file_idx){
	/* File pointer to hold reference to our file */
	FILE * outfile;
	char pname[1024]= "./output_Gurobi/"; 
	char str[1024];
	strcat( pname, "outfile_" );
	sprintf(str, "%s_", gen_str);
	strcat( pname, str );
	sprintf(str, "%d", file_idx);
	strcat( pname, str ); 
	strcat( pname, ".mps" );
	//printf("Generating file %s\n",pname);
	outfile= fopen( pname, "wt" );
	fprintf(outfile, "OBJSENSE    MIN\n");
	fprintf(outfile, "NAME          PROB_%s%d\n", gen_str, file_idx);
	fprintf(outfile, "ROWS\n");
	fprintf(outfile, " E  DEM\n"); // Demand constraint
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile, " L  CAP%d\n",i); // Capacity constraint on node i
	}
	fprintf(outfile, " N  OBJ\n"); // Objective function
	fprintf(outfile, "COLUMNS\n");
	fprintf(outfile, "    MARK0000  'MARKER'                 'INTORG'\n"); // All variables are integer (begin)
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile, "    X%-8d DEM                 1.   CAP%-6d      1.   \n",i,i); // Variable x_i has a coefficient of 1.0 in its capacity constraint and in the demand constraint
		fprintf(outfile, "    X%-8d OBJ       %3.10f\n",i,nodes[i].c);                     // Variable x_i has a coefficient of c_i in the objective  
	}
	fprintf(outfile, "    MARK0001  'MARKER'                 'INTEND'\n"); // All variables are integer (end)
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile, "    Y%-8d CAP%-6d -%-13d\n",i,i,nodes[i].b);  // Variable y_i has a coefficient of -b_i in the ith capacity constraint
		fprintf(outfile, "    Y%-8d OBJ      %3.10f\n",i,nodes[i].f);                     // Variable y_i has a coefficient of f_i in the objective  
	}
	fprintf(outfile, "RHS\n");
	fprintf(outfile, "    B         DEM       %d\n", D); //right-hand side of the demand constraint
	fprintf(outfile, "BOUNDS\n");
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile, " UI BND       X%-8d %13d\n",i,nodes[i].b); // Variable x_i is inferior or equal to bi
	}
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile, " BV BND       Y%-8d\n",i); // Variable y_i is binary
	}
	fprintf(outfile, "ENDATA\n");
	fclose( outfile );
}

//produces the files that will be made publicly available (same format as outputFileKlose)
void outputFileFinal(Node * nodes, int n, int D, char * instance_name, int i){
	/* File pointer to hold reference to our file */
	FILE * outfile;
	char pname[1024]= "./instances/"; 
	char str[1024];
	strcat( pname, "" );
	sprintf(str, "%s_", instance_name);
	strcat( pname, str );
	sprintf(str, "%d", i);
	strcat( pname, str ); 
	//printf("Generating file %s\n",pname);
	outfile= fopen( pname, "wt" );
	fprintf(outfile, "%d\n", n );
	fprintf(outfile, "%d\n", D );
	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%d ",nodes[i].b );
	}
	fprintf(outfile,"\n");

	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%15.16f ", nodes[i].c );
	}
	fprintf(outfile,"\n");

	for (int i = 0; i < n; ++i)
	{
		fprintf(outfile,"%d ",nodes[i].f );
	}
	fclose( outfile );
}


Node* read_file(char * instance_path, int *n, int *D, int *lb, int *ub, double *lc, double *uc, int *lf, int *uf)
{
	FILE *file = fopen(instance_path, "r");
  int curr_int;
  double curr_double;

  // read n
  fscanf(file, "%d", &curr_int);
  *n = curr_int;

  // read D
  fscanf(file, "%d", &curr_int);
  *D = curr_int;


  //initialize nodes
  Node *nodes = (Node*)calloc((*n), sizeof(Node));

  // read b_i, i=1,...,n
  for (int i = 0; i < *n; ++i)
  {
  	fscanf(file, "%d", &curr_int);
  	nodes[i].b = curr_int;
  }

  // read c_i, i=1,...,n
  for (int i = 0; i < *n; ++i)
  {
  	fscanf(file, "%lf", &curr_double);
  	nodes[i].c = curr_double;
  }

  // read f_i, i=1,...,n
  for (int i = 0; i < *n; ++i)
  {
  	fscanf(file, "%d", &curr_int);
  	nodes[i].f = curr_int;
  }

  // complete the initialization of the nodes
  for(int i=0; i<*n; ++i){
    nodes[i].i = i;
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

  // global bounds 
  double minC = DBL_MAX;
	double maxC  = 0.0;
	int minF  = INT_MAX;
	int maxF  = 0.0;
	int minB  = INT_MAX;
	int maxB  = 0.0;
	for(int i=0; i<*n; ++i){
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
      if(nodes[i].b < minB){
        minB = nodes[i].b;
      }
      if(nodes[i].b > maxB){
        maxB = nodes[i].b;
      }
  }

  *lb = minB;
  *ub = maxB;
  *lc = minC;
  *uc = maxC;
  *lf = minF;
  *uf = maxF;

  //printNodes(nodes, *n);
  
  /* Close file */
  fclose(file);
  sort_nodes(*n, nodes);
	return nodes;
}







//======================================================= Wrapper functions ======================================================
double solveSSFCTP_detailed(
	Node *nodes,                 // Suppliers of the problem
	Prob *pr,                    // Problem data
	Res *res,                    // Structure that will contain detailed execution results
	int do_posCostsRed,          // Minimal unit cost c_j centered to 0                                    if true
	int do_sort,                 // Phase 0 executed                                                       if true
	int do_ZLP_ZG,               // Phase 1 executed                                                       if true
	int do_P3,                   // Phases 2 and 3 executed                                                if true
	int do_dom,                  // Phase 6 exectured                                                      if true
	int do_LP,                   // Phase 7 executed                                                       if true
	int do_filter,               // Phase 9 executed                                                       if true
	int do_exact,                // Phases 10 executed                                                     if true
	int computeSol,              // Computes the optimal solution                                          if true, just the optimal value otherwise
	int verbose,                 // Intermediate results printed on standard output                        if true
	int outputSol,               // Prints the optimal solution on standard output                         if true
	int alternative_lambda_init  // Use max_{j in {1,2,...,s}}{c_j} as initial value of lambda in P3Search if true, and min{e_s, max_{j in {1,2,...,n}}{c_j}} otherwise
	){
	//*** Arguments verification ***//
	if(do_posCostsRed != 0 && do_posCostsRed != 1){
		printf("Error: do_posCostsRed must have a binary value (0 or 1)\n"); 
		return 0 ;
	}
	if(do_sort != 0 && do_sort != 1){
		printf("Error: do_sort must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_ZLP_ZG != 0 && do_ZLP_ZG != 1){
		printf("Error: do_ZLP_ZG must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_P3 != 0 && do_P3 != 1){
		printf("Error: do_P3 must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_dom != 0 && do_dom != 1){
		printf("Error: do_dom must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_LP != 0 && do_LP != 1){
		printf("Error: do_LP must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_filter != 0 && do_filter != 1){
		printf("Error: do_filter must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(do_exact != 0 && do_exact != 1){
		printf("Error: do_exact must have a binary value (0 or 1)\n"); 
		return 0;
	}
	if(verbose != 0 && verbose != 1 && verbose != 2){
		printf("Error: verbose must have a value in {0, 1, 2}\n");
		return 0; 
	}

	//*** Variables and structures initialization ***//
	struct timespec start, end0, end1, end5, end6, end7, end9, end;
	P3struct *P3L;
	//res = (Res*)calloc(1, sizeof(Res));
	//res = (Res*) malloc(sizeof(Res));
	int n = pr->n;
	int solvedExaclty = 0; // set to 1 when the optimality of the incumbent solution is proven
	int max_size_L = 10;     //maximum number of knapsacks to solve for the new heuristic

	//*** Positive costs transformation ***//
	clock_gettime(CLOCK_MONOTONIC, &start); // Time is first measured before the positive costs transformation // 
	double fixedObjectiveCost = 0;
	if(do_posCostsRed){
        fixedObjectiveCost = setMinCjToZero(pr, nodes);
        if(verbose==2){
            printf("**************************************************************************************************\n");
            printf("Problem transformation performed. Theta = %f was removed from each node's unit cost. \n", fixedObjectiveCost/pr->D);
            printf("D*theta = %f will be added to the modified problem's optimal objective value at the end of\nthe algorithm to obtain the original  problem's optimal objective value \n", fixedObjectiveCost);
            printf("**************************************************************************************************\n\n");
            printf("%s\n", "Transformed problem");
            printProblem(pr);
            printf("%s\n", "");
        }
    }
	
	//*** Heuristic phase ***//
	// H0) SORTING STEP: 
	// sort the suppliers in non-decreasing order of linearized cost, i.e. e1<=e2<=...<=en
	if(do_sort){
		sort_nodes(n, nodes);
	}
	clock_gettime(CLOCK_MONOTONIC, &end0);
	
	// H1) WEAK LINEAR RELAXATION AND GREEDY PHASE: 
	// calculate the first lower and upper bounds Z_LB and Z_UB given by Z^LP and Z^G
	double incumbentAfter1 = pr->UB;
	if(do_ZLP_ZG){
       	LB_UB_Dantzig_modif(pr, nodes, computeSol);
       	//printf("Classical linear relaxation LB = %f\n", pr->LB);
       	//printf("GREEDY UB                      = %f\n", pr->UB);
       	//printf("GAP                            : %f\n", pr->UB/pr->LB-1);

		if(verbose==2){
			printf("\n%s\n", "Problem after phase 1 (Weak linear relaxation and greedy phase)");
			printProblem(pr);
			printf("%s\n", "");
		}
	    // Optimality check //
		if(pr->LB >= pr->UB){
			solvedExaclty = 1;
			pr->nb_poss_p = 0;
		}
		res->prOpt_ZLP_ZG_1 = solvedExaclty;
		incumbentAfter1 = pr->UB;
	}
    clock_gettime(CLOCK_MONOTONIC, &end1);

	// H2) P3 PHASE: 
	// execute P3Search algorithm to build the set of multiplier values L such that Z^KP(lambda) is known for each lambda in L. Update Z_UB by computing Z^P3_UB
	// if Z_UB <= Z^KP(lambda_0), where lambda_0 is the smallest multiplier value in L, solve P3(0) to restrict future searches to solutions with a partial
	// calculate the lower bounds C^P3_min
	// build the set P of nodes that may improve the incumbent solution while being partial. Compute the bounds LB^P3_xp and UB^P3_xp on each node p in P  
    double incumbentAfter3 = pr->UB;
    if(do_P3 && !solvedExaclty){
    	int kpHeuristicPhase = 0;
	    P3L = P3search(pr, nodes, max_size_L, &kpHeuristicPhase, computeSol, alternative_lambda_init); 
	    incumbentAfter3 = pr->UB;

	    if(pr->LB >= pr->UB){ //if the optimal solution have been found and proven optimal in P3Search, no need to compute CminCmaxLBFilter_LBUpdate, LBp1_all and LB_UBxp1_all
	    	solvedExaclty = 1;
	    	pr->nb_poss_p = 0;
	    }
	    else{ //the optimal solution has not been found and proven optimal in P3Search. Hence, we need to compute CminCmaxLBFilter_LBUpdate, LBp1_all and LB_UBxp1_all
		    CminCmaxLBFilter_LBUpdate(pr, nodes); 
	      if(verbose==2){
	      	printf("\n%s\n", "Execution of P3Search algorithm");
	      	printP3L(P3L, max_size_L);
	      }
	    	LBp1_all(nodes, P3L, pr, n, max_size_L);
	    	LB_UBxp1_all(nodes, P3L, pr, n, max_size_L);
	      CminCmaxLBFilter_LBUpdate(pr, nodes);

	      if(pr->LB >= pr->UB){
		    	solvedExaclty = 1;
		    	pr->nb_poss_p = 0;
		    }
      }

      if(verbose==2){
      	printf("\n%s\n", "Problem after phase 5 (P3 phase and filtering of the partial nodes)");
      	printProblem(pr);
      	printf("%s\n", "");
      }
	    
	    res->kpHeuristicPhase = kpHeuristicPhase;
	    res->prOpt_P3_2_3_4_5 = solvedExaclty;
	    res->cardP_P3_2_3_4_5 = pr->nb_poss_p;
	}
    clock_gettime(CLOCK_MONOTONIC, &end5);

	// H3) DOMINANCE PHASE: 
	// compute LB^Dom_xp for each node
	if(do_dom && !solvedExaclty){
        LBxp3_all_efficient(nodes, n, pr);
        CminCmaxLBFilter_LBUpdate(pr, nodes);
        if(verbose==2){
        	printf("\n%s\n", "Problem after phase 6 (dominance phase)");
        	printProblem(pr);
        	printf("%s\n", "");
        }
        if(pr->LB >= pr->UB){
        	solvedExaclty = 1;
        	pr->nb_poss_p = 0;
        }
        res->prOpt_dom_6 = solvedExaclty;
	    res->cardP_dom_6 = pr->nb_poss_p;
    }
    clock_gettime(CLOCK_MONOTONIC, &end6);

	// H4) STRONG LINEAR RELAXATION PHASE: 
	// compute Z^LP_LB(p), LB^LP_xp and UBLP_xp for each node
	if(do_LP && !solvedExaclty){
        LBp2_LB_UBxp2_all_efficient(nodes, pr, n); 
        CminCmaxLBFilter_LBUpdate(pr, nodes);
        if(verbose==2){
        	printf("%s\n", "Problem after phase 7 (Strong linear relaxation phase)");
        	printProblem(pr);
        	printf("%s\n", "");
        }
        if(pr->LB >= pr->UB){
        	solvedExaclty = 1;
        	pr->nb_poss_p = 0;
        }
        res->prOpt_LP_7 = solvedExaclty;
    	res->cardP_LP_7 = pr->nb_poss_p;
    }
    clock_gettime(CLOCK_MONOTONIC, &end7);


	//*** Filtering phase ***//
	// F) FILTERING PHASE
	// apply dominance filtering and strong linear filtering
	if(do_filter && !solvedExaclty){
        linear_filtering_improved_Q(pr, nodes); 
        //linear_filtering(pr, nodes); 
    }
    clock_gettime(CLOCK_MONOTONIC, &end9);

  //*** Exact phase ***//
	// E) EXACT RESOLUTION 
	// solve the knapsack subproblems P4p, solutions to (P1|p is partial) in \tilde{S}^1 using \tilde{f}_{3,1}
	if(do_exact && !solvedExaclty){
  	if(verbose==2){
  		printf("\n%s\n", "Execution of exact resolution phase\n----------------------------\n");
  	}
  	int kpExactPhase = 0;
  	int totNodesKpExactPhase = 0;
      solveKnapsackFixedPartial_all(nodes, pr, &kpExactPhase, verbose, &totNodesKpExactPhase, computeSol);
      res->kpExactPhase = kpExactPhase;
      res->totNodesKpExactPhase = totNodesKpExactPhase;

      if(verbose==2){
      	printf("\n%s\n", "Problem after phase 10 (Exact resolution)");
      	printProblem(pr);
      	printf("%s\n", "");
      }
      res->prOpt_Exact_10 = 1;
  }
  clock_gettime(CLOCK_MONOTONIC, &end);

  //*** Optimal solution ***//
  if(outputSol){
  	printSolution(pr, nodes);
  }

  //*** Optimal value adjusted if positive costs transformation was performed ***//
  double optimalValue = pr->UB + fixedObjectiveCost;
  if(verbose>=1){
  	printf("Optimal value = %f\n", optimalValue);
  }

  //*** Collect some data ***//
  if(do_exact){ //We can only be sure that the optimal solution has been found if do_exact==1
  	res->isOpt_ZLP_ZG_1 = (incumbentAfter1 == pr->UB);
  	res->isOpt_P3_2_3_4_5 = (incumbentAfter3 == pr->UB);
  }
  res->excessSupply     = getExcessSupply(pr, nodes);
  res->cpu_sort_0       = (end0.tv_sec - start.tv_sec) * 1000000 + (end0.tv_nsec - start.tv_nsec) / 1000; // microseconds;  
  res->cpu_ZLP_ZG_1     = (end1.tv_sec - end0.tv_sec)  * 1000000 + (end1.tv_nsec - end0.tv_nsec)  / 1000;  
 	res->cpu_P3_2_3_4_5   = (end5.tv_sec - end1.tv_sec)  * 1000000 + (end5.tv_nsec - end1.tv_nsec)  / 1000;
 	res->cpu_dom_6        = (end6.tv_sec - end5.tv_sec)  * 1000000 + (end6.tv_nsec - end5.tv_nsec)  / 1000;
 	res->cpu_LP_7         = (end7.tv_sec - end6.tv_sec)  * 1000000 + (end7.tv_nsec - end6.tv_nsec)  / 1000;
 	res->cpu_filtering_9  = (end9.tv_sec - end7.tv_sec)  * 1000000 + (end9.tv_nsec - end7.tv_nsec)  / 1000;
 	res->cpu_exact_10     = (end.tv_sec  - end9.tv_sec)  * 1000000 + (end.tv_nsec  - end9.tv_nsec)  / 1000;
 	res->cpu_tot          = (end.tv_sec - start.tv_sec)  * 1000000 + (end.tv_nsec  - start.tv_nsec) / 1000;
 	res->opt_val          = optimalValue;
 	if(verbose>=1){
 		printf("Total CPU time = %f ms\n", (double) ((double)res->cpu_tot/1000));
 	}

 	return optimalValue;
}

double solveSSFCTP(
	Node *nodes,   // Suppliers of the problem
	Prob *pr,      // Problem data
	Res *res       // Structure that will contain detailed execution results
	){
	return solveSSFCTP_detailed(nodes, pr, res, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0);
}

//return the optimal value of the SSFCTP(b,c,f,D,n) and its optimal solution (x,y)
//this function supports negative values for both variable costs c_i and 
//fixed costs f_j
double kaa_solve(int * b, double * c, double * f, int D, int n, int * x, int * y){

	//Deals with negative fixed costs f_j
	double sum_negative_f = 0;
	int * is_f_neg = (int*)calloc(n, sizeof(int));
	double * f_non_neg = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; ++i){
		f_non_neg[i] = f[i];
		if(f[i] < 0){
			sum_negative_f += f[i];
			is_f_neg[i] = 1;
			f_non_neg[i] = 0;
		}
	}

	Node *nodes = get_nodes(b, c, f_non_neg, n);
	//printNodes(nodes, n);
	Prob *pr = get_prob(D, b, c, f_non_neg, n);
	Res *res = create_res();
	double Z = solveSSFCTP(nodes, pr, res);
	for (int i = 0; i < n; ++i){
		int idx = nodes[i].i_init;
		x[idx] = pr->incumbent[i];
		y[idx] = (pr->incumbent[i] || is_f_neg[idx]); //y_i=1 <=> (positive usage x_i or negative fixed cost f_i)
	}

	return Z + sum_negative_f;
}



