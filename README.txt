Optimal transport with proximal splitting - MATLAB CODE

G. PeyrŽ, N. Papadakis, E. Oudet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This archive contains the Douglas-Rachford (DR) and Primal-Dual (PD) solvers applied to the Benamou-Brenier (BB) problem
discretized on a staggered grid.

They can be tested with:

test_bb_dr
test_bb_pd

%%%%%%%%%%%%%%%%%%%%

the principal options are the following:

1) chose your test case with 
test = 'gaussian';
(you can create your own scenario by defining f0 and f1 the initial and final densities)

2) chose the dimension of the problem:
N=32; P=32; Q=32;  (N and P are the discrete spatial dimensions and Q is the temporal discretization)


3) Parameterization of the solver:

DR: 
mu = 1.98; % should be in ]0,2[
gamma = 1./230.; % should be >0
niter = 1000;


PD:
sigma=85;
niter = 1000;   (increase the maximum number of iteration to have better results)

4) Generalized cost functions:

Minimize \sum_k w_k f_k^\alpha |v_k|^2
alpha= 1; % should be in [0;1];

%alpha=1 computes the L2-Wasserstein distance, 0 is for the H^-1 one and intermediate values gives interpolations between the norms

obstacle=zeros(N,P,Q);
%define the points of the 3D volume where the mass can not pass 
For instance, setting 
obstacle(N/2,P/2,:)=1;
will create an obstacle in the middle of the spatial domain.

%%%%%%%%%%%%%%%%%%%%%%
Exemples:

test = 'gaussian';
N=32; P=32; Q=32;  
niter = 200;


test = 'obstacle';
niter = 2000;