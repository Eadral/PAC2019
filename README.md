# PAC2019优化组决赛  手算队201



优化内容：

1. 算法
2. 编译选项
3. 杂项



> This readme file contains compilation and run scripts for the sample codes.
>                                                             Bin 10/22/2015
> pcg_sfc_quiz: 
> compile: compile_pcg.sh
> run    : run_pcg.sh
> 
> Notices:
> 1) Residual threshold is adjusted so that double-precision floating-point might 
> be necessary for convergence.
> 2) MPI rank number is adjusted to 24 so that no core is oversubscribed on HSW/
>  BDW/SKL/CLX processor.
> 3) The result should be reproducible. Please submit all changes, including codes, 
> building and running scripts.
> 4) Can't modified the file: pcg.90 , solver_pcg_mod.f90 , time_management.f90, pop_in
> 4) Algorithm optimization is strongly encouraged. Some references on Communication-
> avoiding Krylov subspace methods are given below, as well as the P-CSI (Chebyshev 
> Iteration) paper.
> Some infrastructural functions (hints) are already provided, such as create_boundary_capcg 
> and update_matpow_halo for thicker boundary update in CA-PCG solver, global_sum_1d_dbl 
> for one-dimensional double array, etc.
> 
> References:
> [1] M. Hoemmen, Communication-avoiding Krylov subspace methods, PhD thesis,
> University of California, Berkeley, 2010.
> [2] Erin Carson, Communication-Avoiding Krylov Subspace Methods in Theory and
> Practice, PhD thesis, University of California, Berkeley, 2015.
> [3] Erin Carson and James Demmel, A Residual Replacement Strategy for Improving
> the Maximum Attainable Accuracy of s-step Krylov Subspace Methods.
> SIAM J. Matrix Anal. Appl., 35(1), 2014, pp. 22-43.
> [4] J. Demmel, M. Hoemmen, M. Mohiyuddin, and K. Yelick, Avoiding communication
> in computing Krylov subspaces, Tech. Report UCB/EECS-2007-123, EECS Dept., U.C.
> Berkeley, Oct 2007.
> [5] Yong Hu, Xiaomeng Huang, etal. Improving the Scalability of the Ocean Barotropic 
> Solver in the Community Earth System Model. SC 2015.

