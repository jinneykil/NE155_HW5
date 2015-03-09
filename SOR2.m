function [x, k]=SOR2(a, b, x0, tol, w)
% This function uses the iterative method: Successive Over Relaxation (SOR)
% INPUT
% a: matrix A in Ax=b
% b: matrix b in Ax=b
% x0: initial x^(0)
% tol: absolute tolerance
% w: scalar w
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
            sum1(j)=-(w*a(i,j)/a(i,i))*x(j);
        end
    end
    xn=x;
    x(i)=sum(sum1)+(w*b(i)/a(i,i))+((1-w)*(x(i)));
    sum1=0;
end
k=k+1;
    
while tol<(norm(xn-x)/norm(xn))
    for i=1:n
        for j=1:n
            if j==i
                sum1(j)=0;
            else
                sum1(j)=-(w*a(i,j)/a(i,i))*x(j);
            end
        end
        xn=x;
        x(i)=sum(sum1)+(w*b(i)/a(i,i))+((1-w)*(x(i)));
        sum1=0;
    end
    k=k+1;
end