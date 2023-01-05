#include "util.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "ssfctp.h" //TODO: à enlever (juste nécessaire pour void printNodesPart() actuellement.)
                    //      les fonctions propres au ssfctp (build_heap_LR...) devront être déplacées

//============================================= Nb manipulations ====================================================
int ceil_log2(int x){
    return ceil(log(x) * M_LOG2E);
}

//======================================== Random number generation =================================================

//return a random double uniformly in (0,1)
double randU01()
{
	return (double)rand() / (double)RAND_MAX ;
}

//return a random integer uniformly in {min, min+1,... , max}
int randInt(int min, int max){
	return (int) (min + (max-min)*randU01());
}

//return a random double uniformly in (min,max)
double randDouble(double min, double max){
	return (min + (max-min)*randU01());
}

//TODO: index_comparator et partial_node_comparator font la même chose actuellement, sauf pour des types d'entrée diff.
// return -1 if n1<n2 according to the article's total order (prop. 11.2), i.e.
// (c_1<c_2) or ((c1=c2) and (b1<b2)) or ((c1=c2) and (b1=b2) and (i1<i2)). Return +1 otherwise. 
// If index_comparator(*n1, *n2)>0, node 1 is excluded from any solution using node 2 as a partial.
// n1 and n2 are pointers to pointers to nodes
static int index_comparator(const void *n1, const void *n2)
{
	const Node *p1 = *(Node **)n1;
	const Node *p2 = *(Node **)n2;
	if ((p1->c<p2->c)||((p1->c==p2->c)&&(p1->b<p2->b))||((p1->c==p2->c)&&(p1->b==p2->b)&&(p1->i_init<p2->i_init)))  //n1<n2. i.e. 1 partial => 2 unused
		return -1;
	return +1;
}

// return  0 if n1=n2 (i.e. same node passed twice), 
// return -1 if n1<n2 (i.e. 1 partial => 2 unused), 
// return +1 if n1>n2 (i.e. 2 partial => 1 unused)
int partial_node_comparator(Node *p1, Node *p2)
{
	if(p1->i_init==p2->i_init) //n1=n2.
		return 0;
	if ((p1->c<p2->c)||((p1->c==p2->c)&&(p1->b<p2->b))||((p1->c==p2->c)&&(p1->b==p2->b)&&(p1->i_init<p2->i_init)))  //n1<n2. 
		return -1;
	return +1;
}

//=========================== Heap functions (to be used on an array of pointers) ===================================

//insertion of a new pointer to a node in maxHeap nodesLR
int heap_insert(Node** nodesLR, Node* ptrNode, int* nbElementsLR)
{
	int index;
	int parent;

	/* Add to the bottom of the heap and start from there */
	index = *nbElementsLR;
	*nbElementsLR += 1;

	/* Percolate the value up to the top of the heap */
	while (index > 0) {

		/* The parent index is found by halving the node index */
		parent = (index - 1) / 2;

		/* Compare the node with its parent */
		if (partial_node_comparator(nodesLR[parent], ptrNode) > 0) {

			/* Ordered correctly (the parent is bigger than ptrNode according to prop. 11.2 order) - insertion is complete */
			break;
		} else {

			/* Need to swap this node with its parent */
			nodesLR[index] = nodesLR[parent];

			/* Advance up to the parent */
			index = parent;
		}
	}

	/* Save the new value in the final location */
	nodesLR[index] = ptrNode;
	return 1;
}

Node* heap_pop(Node** nodesLR, int* nbElementsLR)
{
	Node* result;
	Node* new_value;
	int index, next_index, child1, child2;

	/* Empty heap? */
	if (*nbElementsLR == 0) {
		return NULL;
	}

	/* Take the value from the top of the heap */
	result = nodesLR[0];

	/* Remove the last value from the heap; we will percolate this down from the top. */
	new_value = nodesLR[*nbElementsLR-1];
	*nbElementsLR-=1;

	/* Percolate the new top value down */
	index = 0;

	for (;;) {

		/* Calculate the array indexes of the children of this node */
		child1 = index * 2 + 1;
		child2 = index * 2 + 2;

		if (child1 < *nbElementsLR && partial_node_comparator(nodesLR[child1], new_value) > 0) {

			/* Left child is bigger than the node.  We need to swap
			 * with one of the children, whichever is less. */
			if (child2 < *nbElementsLR && partial_node_comparator(nodesLR[child2], nodesLR[child1]) > 0) {
				next_index = child2;
			} 
			else {
				next_index = child1;
			}

		} 
		else if (child2 < *nbElementsLR && partial_node_comparator(nodesLR[child2], new_value) > 0) {

			/* Right child is bigger than the node.  Swap with the right child. */
			next_index = child2;

		} 
		else {

			/* Node is bigger than both its children. The heap
			 * condition is satisfied.  * We can stop percolating
			 * down. */
			nodesLR[index] = new_value;
			break;
		}

		/* Swap the current node with the least of the child nodes. */
		nodesLR[index] = nodesLR[next_index];

		/* Advance to the child we chose */
		index = next_index;
	}

	return result;
}

void printHeap(Node** nodesLR, int* nbElementsLR){
	for (int i = 0; i < *nbElementsLR; ++i)
	{
		printf("%f\n",nodesLR[i]->c);
	}
}

//======================================== Build functions for LBp2 =================================================

// Initial build of the max-heap nodesLR, by including the nodes used in P3's regular linear relaxation.
// The objective value ZLR, the split node index s and the number of elements nbElementsLR in nodesLR
// are stored in the correspondant variables (passed as pointers).
void build_heap_LR(int n, int D, double *ZLR, int *s, Node *nodes, Node **nodesLR, int *nbElementsLR){
	for (int i = 0; i < n; ++i)
	{
		heap_insert(nodesLR, &nodes[i], nbElementsLR);
		if(D > nodes[i].b){
			nodes[i].x = nodes[i].b;
			D -= nodes[i].x;
			*ZLR += nodes[i].u;
		}
		else{
			nodes[i].x = D;
			D = 0; // <=> D -= nodes[i].x;
			*ZLR += nodes[i].x * nodes[i].e;
			*s = i; //split node's index
			break;
		}
	}
}

// Builds the linked list nodesLinked containing the problem's nodes, sorted in non-decreasing order of linearized cost e_j (it is assumed that
// the array nodes respects this order). firstNodeLL, splitNodeLL and lastNodeLL are pointers to the corresponding NodeLinkedList.
void build_linkedList_nodes(NodeLinkedList **firstNodeLL, NodeLinkedList **splitNodeLL, NodeLinkedList **lastNodeLL, Node *nodes, int s, int n){
	NodeLinkedList *currentNodeLL;
	(*firstNodeLL) = (NodeLinkedList*)malloc(sizeof(NodeLinkedList)); // first nodeLL
	(*firstNodeLL)->prev = NULL;                                      //  - no previous nodeLL
	(*firstNodeLL)->ptrNode = &nodes[0];                              //  - data: &nodes[0]
	(*splitNodeLL) = (*firstNodeLL); //split nodeLL (points to nodes[0] initially)

	currentNodeLL = (*firstNodeLL);
	for (int i = 1; i < n; ++i) //intermediate nodes
	{
		NodeLinkedList *newNodeLL = (NodeLinkedList*)malloc(sizeof(NodeLinkedList)); //new nodeLL
		currentNodeLL->next= newNodeLL;                                              //  - connect current                
		newNodeLL->prev = currentNodeLL;                                             //    and new nodesLL
		newNodeLL->ptrNode = &nodes[i];                                              //  - data: &nodes[i]
		currentNodeLL=newNodeLL;
		if(currentNodeLL->ptrNode->i == s){//split nodeLL: now points to the real split node
			(*splitNodeLL) = currentNodeLL;   //(if s!=0; otherwise it will continue to point to 0)
		}
	}

	(*lastNodeLL) = currentNodeLL; // last nodeLL
	(*lastNodeLL)->next = NULL;    // -no next nodeLL
}

//Builds an array list of pointers to nodes that may be partial in an optimal solution to P3, sorted in 
//			  increasing order (using total order from the article)
void build_arraylist_partials(int n, Node * nodes, Node **nodesPart, int *nb_poss_part){
	*nb_poss_part = 0;
	for (int i = 0; i < n; ++i)
	{
		if(nodes[i].poss_p){
			nodesPart[*nb_poss_part] = &nodes[i];
			*nb_poss_part += 1;
		}
		
	}
	qsort(nodesPart, *nb_poss_part, sizeof(nodesPart[0]), index_comparator);
}

//Free the memory space used by the linked list with first and last nodes firstNodeLL, lastNodeLL
void util_free_LL(NodeLinkedList* firstNodeLL, NodeLinkedList* lastNodeLL){
	NodeLinkedList* currentNodeLL;
	if(firstNodeLL != NULL){//free the linkedlist
		currentNodeLL = firstNodeLL; 
		while(currentNodeLL != NULL){ 
			NodeLinkedList *temp = currentNodeLL->next;
			free(currentNodeLL);
			currentNodeLL = temp;
		}
	}
	else if(lastNodeLL != NULL){
		currentNodeLL = lastNodeLL;
		while(currentNodeLL != NULL){
			NodeLinkedList *temp = currentNodeLL->prev;
			free(currentNodeLL);
			currentNodeLL = temp;
		}
	}
	firstNodeLL = NULL;
	lastNodeLL = NULL;
}



//======================================================= Output functions =======================================================
void printNodesPart(Node **nodesPart, int nb_poss_part){
	for (int i = 0; i < nb_poss_part; ++i)
	{
		printNodeDetailed(nodesPart[i],0);
	}
}

//============================================ Detailed results container operations =============================================
GlobalRes* createGlobalRes(){
	GlobalRes* gRes = (GlobalRes*)calloc(1, sizeof(GlobalRes));
	return gRes;
}

Res* create_res(){
	Res* res = (Res*)calloc(1, sizeof(Res));
	return res;
}

void cleanMemory(Res* res, P3struct* P3L, Prob* pr, Node* nodes){
	if(res->prOpt_ZLP_ZG_1 == 0){ //this means that P3Search algorithm has been performed, hence P3L has been initialized
		free(P3L);
	}
	free(pr);
	free(nodes);
}

void updateGlobalRes(GlobalRes *gRes, Res * res){
	gRes->cpu_sort_0           += res->cpu_sort_0;
	gRes->cpu_ZLP_ZG_1         += res->cpu_ZLP_ZG_1;
	gRes->cpu_P3_2_3_4_5       += res->cpu_P3_2_3_4_5;
	gRes->cpu_dom_6            += res->cpu_dom_6;
	gRes->cpu_LP_7             += res->cpu_LP_7;
	gRes->cpu_filtering_9      += res->cpu_filtering_9;
	gRes->cpu_exact_10         += res->cpu_exact_10;
	gRes->cpu_tot              += res->cpu_tot;
	gRes->prOpt_ZLP_ZG_1       += res->prOpt_ZLP_ZG_1;
	gRes->isOpt_ZLP_ZG_1       += res->isOpt_ZLP_ZG_1;
	gRes->prOpt_P3_2_3_4_5     += res->prOpt_P3_2_3_4_5;
	gRes->isOpt_P3_2_3_4_5     += res->isOpt_P3_2_3_4_5;
	gRes->prOpt_dom_6          += res->prOpt_dom_6;
	gRes->prOpt_LP_7           += res->prOpt_LP_7;
	gRes->prOpt_Exact_10       += res->prOpt_Exact_10;
	gRes->kpHeuristicPhase     += res->kpHeuristicPhase;
	gRes->kpExactPhase         += res->kpExactPhase;
	gRes->totNodesKpExactPhase += res->totNodesKpExactPhase;
	gRes->cardP_P3_2_3_4_5     += res->cardP_P3_2_3_4_5;
	gRes->cardP_dom_6          += res->cardP_dom_6;
	gRes->cardP_LP_7           += res->cardP_LP_7;
	gRes->excessSupply         += res->excessSupply;  
    gRes->numberPartials       += (res->excessSupply >= 1);
	gRes->nRep                 += 1;
	gRes->sum_opt_val          += res->opt_val;
}

void printGlobalRes(GlobalRes *gRes, int computeSol){
	printf("\n========================================= GLOBAL STATISTICS =========================================\n");
	printf("Number of instances        : %d\n", gRes->nRep);
	printf("Mean cpu time per instance : %f%s\n", (double)gRes->cpu_tot/1000/gRes->nRep, " msec");
	printf("Sum of optimal values      : %f\n", gRes->sum_opt_val);
	printf("\n_____________________________________________ CPU TIME ______________________________________________\n");
	printf("Cpu time in phase  H1 :  %f%s\n", (double)gRes->cpu_sort_0+(double)gRes->cpu_ZLP_ZG_1 /1000, " msec");
	printf("Cpu time in phase  H2 :  %f%s\n", (double)gRes->cpu_P3_2_3_4_5 /1000, " msec");
	printf("Cpu time in phase  H3 :  %f%s\n", (double)gRes->cpu_dom_6      /1000, " msec");
	printf("Cpu time in phase  H4 :  %f%s\n", (double)gRes->cpu_LP_7       /1000, " msec");
	printf("Cpu time in phase  F  :  %f%s\n", (double)gRes->cpu_filtering_9/1000, " msec");
	printf("Cpu time in phase  E  :  %f%s\n", (double)gRes->cpu_exact_10   /1000, " msec");
	printf("Total cpu time        :  %f%s\n", (double)gRes->cpu_tot        /1000, " msec");
	printf("\n______________________________________ RESOLUTION STATISTICS ________________________________________\n");
	printf("Instances solved exactly during phase  H1  : %d\n", gRes->prOpt_ZLP_ZG_1);
	printf("Instances solved exactly during phases H2  : %d\n", gRes->prOpt_P3_2_3_4_5);
	//printf("Instances solved exactly during phase H3 : %d\n", gRes->prOpt_dom_6); //(cannot happen, since no new solution is identified during this phase)
	printf("Instances solved exactly during phase  H4  : %d\n", gRes->prOpt_LP_7);
	printf("Instances solved exactly during phase  E   : %d\n", gRes->prOpt_Exact_10);
	printf("Total                                      : %d\n", gRes->prOpt_ZLP_ZG_1 + gRes->prOpt_P3_2_3_4_5 + gRes->prOpt_LP_7 + gRes->prOpt_Exact_10);
	if(computeSol){
		printf("Number of instances for which the optimal solution contains a partial node : %d \n",gRes->numberPartials);
		if(gRes->numberPartials>=1){
			printf("Mean number of excess items on the partial nodes                           : %f \n",(double)gRes->excessSupply/(double)gRes->numberPartials);
		}
		else{
			printf("Mean number of excess items on the partial nodes                           : %f \n",0);
		}
	}
	printf("\n_______________________________________ P SUBSET STATISTICS _________________________________________\n");
	printf("Mean cardinality of set P after H2 : %f\n", (double)gRes->cardP_P3_2_3_4_5/gRes->nRep);
	printf("Mean cardinality of set P after H3 : %f\n", (double)gRes->cardP_dom_6/gRes->nRep);
	printf("Mean cardinality of set P after H4 : %f\n", (double)gRes->cardP_LP_7/gRes->nRep);
	printf("\n_________________________________ KNAPSACK SUBPROBLEMS STATISTICS ___________________________________\n");
	printf("Mean number of knapsack problems solved during the heuristic phases : %f\n", (double)gRes->kpHeuristicPhase/gRes->nRep);
	printf("Mean number of knapsack problems solved during the exact phase      : %f\n", (double)gRes->kpExactPhase/gRes->nRep);
	if(gRes->prOpt_Exact_10 == 0){
		printf("Mean number of knapsack problems solved during the exact phase when it is executed : %f\n", 0);
		printf("Mean number of items per knapsack problem after the filtering procedure            : %f\n", 0);
	}
	else{
		printf("Mean number of knapsack problems solved during the exact phase when it is executed : %f\n", (double)gRes->kpExactPhase/gRes->prOpt_Exact_10);
		printf("Mean number of items per knapsack problem after the filtering procedure            : %f\n", (double)gRes->totNodesKpExactPhase/gRes->kpExactPhase);
	}
	printf("Mean number of knapsack problems solved during the complete algorithm              : %f\n", (double)(gRes->kpHeuristicPhase+gRes->kpExactPhase)/gRes->nRep);
	printf("=====================================================================================================\n");
}

void advanceSeed(int * seed){
	int a = 31792125;
    int m = 268435399;
    *seed = (a*(*seed))%m; //LCG to generate subsequent seeds
    *seed += m*((*seed)<0);
}

int initSeed(int seed){
	if(seed != 0){
		return seed;
	}
	else{
		time_t tt;
	    time( &tt );
	    return (int) tt;
	}
}
