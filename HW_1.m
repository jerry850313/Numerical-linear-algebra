clc;
i99 = ones(1,99);
n = -1*i99;
I = eye(100);
A = 2*I+diag(n,1)+diag(n,-1);
b = ones(100,1);
x0 = ones(100,1);
for i= 2 : 99
    b(i) = 0;
end
for i= 2 : 100
   x0(i) = 0;
end
TOL=1e-10;
N0=1e5;
[N error]=Gauss_Seidel_iteration(A,b,x0,TOL,N0);
N
error
semilogy(N,error);
grid on;
