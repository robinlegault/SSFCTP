# The Knapsack Transformation Algorithm (KTA) for solving the single-sink fixed-charge transportation problem (SSFCTP)

The KTA is described in "A Novel Reformulation for the Single-Sink Fixed-Charge Transportation Problem" by R.Legault, JF.Côté, B.Gendron.

### ARTICLE ABSTRACT

The single-sink fixed-charge transportation problem is known to have many applications in the area of Manufacturing and Transportation as well to be an important subproblem of the fixed-charge transportation problem. However, even the best solutions of the literature do not succeed in fully exploiting the structure of this problem, to the point of being surpassed by modern general-purpose mixed-integer programming solvers for large instances.

We introduce a novel reformulation of the problem and study its theoretical properties, which lead to a range of new upper and lower bounds as well as specific dominance relations, linear relaxations and filtering procedures. The resulting algorithm includes a heuristic phase and an exact phase, the main step of which is to solve a very small number of knapsack subproblems.

Computational experiments are presented for existing and new types of data instances. These tests indicate that the new algorithm systematically reduces the resolution time of the state-of-the-art exact methods by several orders of magnitude. 


### DIRECTORY STRUCTURE

`main.c` : example programs showing how to generate instances of the SSFCTP and how to solve them using the KTA 

`Makefile` : simple Makefile that can be used to execute the main file given as an example

`combo.c` : slightly modified version of the combo algorithm code presented in "Dynamic programming and tight bounds
for the 0-1 knapsack problem". The original code is made available by the authors: http://hjemmesider.diku.dk/~pisinger/codes.html.
The combo algorithm is used within the KTA to solve the knapsack subproblems.

`combo.h` : header file associated with `combo.c`

`generator.c` : implements a range of classical and new instance generation methods for the SSFCTP

`generator.h` : header file associated with `generator.c`

`ssfctp.c` : principal file. Includes the functions `solveSSFCTP` and `solveSSFCTP_detailed` that can be used to solve instances
of the SSFCTP using the KTA

`util.c` : defines basic operations on data structures, random number generation method, etc.

`util.h` : header file associated with `util.c`. Defines the custom data structures that are used in the rest of the code
