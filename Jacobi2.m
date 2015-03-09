function [x,k]=Jacobi2(a, b, x0, tol)
% This function uses the iterative method: Jacobi
% INPUT
% a: matrix A in Ax=b
% b: matrix b in Ax=b
% x0: initial x^(0)
% tol: absolute tolerance
% OUTPUT
% x: solution
% k: number of iterations to meet tolerance (tol)

n=length(b);    % number of unknowns
sum1=zeros(n,1);
x=x0;
xn=x0;
k=0;    % initializes iteration counter

for i=1:n
    for j=1:n
        if j==i
            sum1(j)=0;
        else
            sum1(j)=-(a(i,j)/a(i,i))*x0(j);
        end
    end
    xn=x;
    x(i)=sum(sum1)+(b(i)/a(i,i));
    sum1=0;
end
x0=x;
k=k+1;

while tol<(norm(xn-x)/norm(xn))
for i=1:n
    for j=1:n
        if j==i
            sum1(j)=0;
        else
            sum1(j)=-(a(i,j)/a(i,i))*x0(j);
        end
    end
    xn(i)=x(i);
    x(i)=sum(sum1)+(b(i)/a(i,i));
    sum1=0;
end
x0=x;
k=k+1;
end
