%% Question 4 Part A
clear all;
close all;
clc;

A=[4 -1 0 0
   -1 4 -1 0
   0 -1 4 -1
   0 0 -1 4 ];
b=[100
    100
    100
    100];
x0=[0
    0
    0
    0];
tol1=10e-6;
tol2=10e-8;
w=1.1;

% Tolerance = 10e-6
[x1,k01]=GaussSeidel2(A,b,x0, tol1);
[y1, k02]=Jacobi2(A, b, x0, tol1);
[z1,k03]=SOR2(A, b, x0, tol1, w);

[x2,k11]=GaussSeidel(A,b,x0, tol1);
[y2, k12]=Jacobi(A, b, x0, tol1);
[z2,k13]=SOR(A, b, x0, tol1, w);

% Tolerance = 10e-8
[x1,k1]=GaussSeidel2(A,b,x0, tol2);
[y1, k2]=Jacobi2(A, b, x0, tol2);
[z1,k3]=SOR2(A, b, x0, tol2, w);

[x2,k21]=GaussSeidel(A,b,x0, tol2);
[y2, k22]=Jacobi(A, b, x0, tol2);
[z2,k23]=SOR(A, b, x0, tol2, w);

% Number of iterations with tol = 10e-6
    % Row 1: relative error
    % Row 2: absolute error
k_tol1=[k01 k02 k03; k11 k12 k13]

% Number of iterations with tol = 10e-8
    % Row 1: relative error
    % Row 2: absolute error
k_tol2=[k1 k2 k3; k21 k22 k23]

% For each method, how does the absolute error with tol=10e-6 compare to t
% he relative error?
    % The absolute error is greater than the relative error for all
    % methods. 
    % GaussSeidel: ab_error-rel_error = 2
    % Jacobi: ab_error-rel_error = 5
    % SOR: ab_error-rel_error = 3
% Which method required the fewest iterations?
    % SOR required the fewest number of iterations
% What do you observe about reaching a tighter convergence tolerance?
    % Reaching a tighter convergence requires more iterations

%% Question 4 Part B
clear all;
close all;
clc;

A=[4 -1 0 0
   -1 4 -1 0
   0 -1 4 -1
   0 0 -1 4 ];
b=[100
    100
    100
    100];
x0=[0
    0
    0
    0];
[z,w_opt]=SOR3(A, b, x0);
w_opt


