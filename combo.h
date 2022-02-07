
/* ======================================================================
             combo.h,    S.Martello, D.Pisinger, P.Toth     feb 1997
   ====================================================================== */

/* This is the header file of the COMBO algorithm described in
 * S.Martello, D.Pisinger, P.Toth: "Dynamic Programming and Strong
 * Bounds for the 0-1 Knapsack Problem".
 *
 *
 /* Modifications by Robin Legault and Jean-François Côté, October 2020 
 */


#ifndef COMBO
#define COMBO
#include <inttypes.h>

#ifdef __cplusplus
  extern "C" {
#endif


//return the maximum profit and x is an array indicating if the item is in or not (0 or 1)
double combo_solve(double * p, double * w, int * x, int64_t W, double lb, double ub, int n);


#ifdef __cplusplus
  }
#endif

#endif
