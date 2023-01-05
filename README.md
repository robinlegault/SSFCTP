# The Knapsack Transformation Algorithm (KTA) for solving the single-sink fixed-charge transportation problem (SSFCTP)

The KTA is described in "A Novel Reformulation for the Single-Sink Fixed-Charge Transportation Problem" by R.Legault, JF.Côté and B.Gendron. Preprint available: https://www.cirrelt.ca/documentstravail/cirrelt-2022-05.pdf 

### ARTICLE ABSTRACT

The single-sink fixed-charge transportation problem is known to have many applications in the area of manufacturing and transportation as well as being an important subproblem of the fixed-charge transportation problem. However, even the best algorithms from the literature do not fully leverage the structure of this problem, to the point of being surpassed by modern general-purpose mixed-integer programming solvers for large instances. We introduce a novel reformulation of the problem and study its theoretical properties. This reformulation leads to a range of new upper and lower bounds, dominance relations, linear relaxations, and filtering procedures. The resulting algorithm includes a heuristic phase and an exact phase, the main step of which is to solve a very small number of knapsack subproblems. Computational experiments are presented for existing and new types of instances. These tests indicate that the new algorithm systematically reduces the resolution time of the state-of-the-art exact methods by several orders of magnitude.


### DIRECTORY STRUCTURE

`main.c` : simple example showing how to solve an instance of the SSFCTP using KTA

`main_g2.c` : example showing how to solve multiple instances of the SSFCTP using KTA to collect statistics (instances of Group 2)

`Makefile` : simple Makefile that can be used to execute the main file given as an example

`combo.c` : slightly modified version of the combo algorithm presented in "Dynamic programming and tight bounds
for the 0-1 knapsack problem". The original code is made available by the authors: http://hjemmesider.diku.dk/~pisinger/codes.html.
The combo algorithm is used within KTA to solve knapsack subproblems.

`combo.h` : header file associated with `combo.c`

`generator.c` : implements a range of classical and new instance generation methods for the SSFCTP

`generator.h` : header file associated with `generator.c`

`ssfctp.c` : principal file. Includes the functions `solveSSFCTP` and `solveSSFCTP_detailed` that can be used to solve instances
of the SSFCTP using KTA

`util.c` : defines basic operations on data structures, random number generation method, etc.

`util.h` : header file associated with `util.c`. Defines the custom data structures that are used in the rest of the code

`instances` : folder containing the instances used for the experiments of Section 7.2

`results` : folder containing the detailed results of the experiments presented in Section 7.2
