Matlab codes for solving saddle point equations under RS ansatz.

- CS_L1_minimization.m: RS-saddle point equations for L1 minimization in compressed sensing(CS)
- CS_L1_minimization_success.m: Stability of the success solution for L1 minimization in CS

	- csolve.m: quasi-Newton method. Returns rc=0 if the routine converges.
	  Use as "csolve(FUN,x,gradfun,crit,itmax,varargin)"
		* FUN: function to be solved
		* x: initial value
		* gradfun: gradient of the FUN
		* crit: error torelance
		* itmax: maximum value of the iteration
		* varargin: parameters required to solving FUN
	- chih_solve.m: function to be solved
	- chih_jacobi.m: gradient of the function
	- See for details https://iopscience.iop.org/article/10.1088/1742-5468/2009/09/L09003
	  and its erratum https://iopscience.iop.org/article/10.1088/1742-5468/2012/07/E07001

