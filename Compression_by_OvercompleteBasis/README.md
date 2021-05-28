Matlab codes for solving saddle point equations of the compression problem by overcomplete basis.
The codes follow the analysis shown in Nakanishi-Ohno et al, JSTAT (2016)
http://dx.doi.org/10.1088/1742-5468/2016/06/063302

- Complexity.m: RS-saddle point equations. The output is the distortion and complexity as a function of the weight parameter mu.

	- csolve.m: quasi-Newton method. Returns rc=0 if the routine converges.
	  Use as "csolve(FUN,x,gradfun,crit,itmax,varargin)"
		* FUN: function to be solved
		* x: initial value
		* gradfun: gradient of the FUN
		* crit: error torelance
		* itmax: maximum value of the iteration
		* varargin: parameters required to solving FUN
	- Rh_solve.m: function for a conjugate variable appears in saddle point equation
	- Rh_jacobi.m: gradient of the function Rh_solve

