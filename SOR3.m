function [x, w_opt]=SOR3(a, b, x0)
% This function uses the iterative method: Successive Over Relaxation (SOR)
% INPUT
% a: matrix A in Ax=b
% b: matrix b in Ax=b
% x0: initial x^(0)
% OUTPUT
% x: solution
% w_opt: optimal w for the given problem

% Procedure:
% SOR performed for different values of w between 1 and 2 (1<w_opt<2)
% For each index, w value, error, and x vector stored in array
% w_opt is the w value of the minimum error:
    % Minimum value and index of min value found for the error array
    % w_opt = w value with same index as minimum error value

n=length(b);    % number of unknowns
sum1=zeros(n,1);
x=x0;
xn=x0;
k=1;    % initializes iteration counter
error=[];
w=[];

for m=1:0.01:2 % different values for w
    for i=1:n
        for j=1:n
            if j==i
                sum1(j)=0;
            else
                sum1(j)=-(m*a(i,j)/a(i,i))*x(j);
            end
        end
        xn=x;
        x(i,k)=sum(sum1)+(m*b(i)/a(i,i))+((1-m)*(x(i)));
        sum1=0;
    end
    w(k)=m; % array stores value of w used for given index, k
    error(k)=norm(xn-x); % calculates the error for each w
    k=k+1;
end
[MIN_er,I]=min(error(:)); % MIN_er: minimum error, I: index of minimum error
w_opt=w(I); % finds w with same index of min error


