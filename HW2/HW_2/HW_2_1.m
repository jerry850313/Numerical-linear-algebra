function HW_2_1
clc;clear all;
function  [error] = Gauss_Seidel_iteration(A,b,x0)
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
Tg = -(D+L)\U;
Cg = (D+L)\b;
x = Tg*x0+Cg;
t=1;
TOL=1e-10;
error(t)=norm(x-x0);
T=norm(x-x0)/norm(x0);
x0=x;
    while T>=TOL
        x = Tg*x0+Cg;
        error(t+1)=norm(x-x0);
        t=t+1;
        x0=x;
        if t>(1e5-1)
            return
        end
    end
end
function  [error] = jacobi_iteration(A,b,x0)
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
B = D\(L+U);
f = D\b;
x = B*x0+f;
t=1;
TOL=1e-10;
error(t)=norm(x-x0);
T=norm(x-x0)/norm(x0);
x0=x;
    while T>=TOL
        x = B*x0+f;
        error(t+1)=norm(x-x0);
        someone=error(t+1);
        t=t+1;
        x0=x;
        if t>(1e5-1)
            return
        end
    end
end
i99 = ones(1,99);
n = -1*i99;
I = eye(100);
A = 2*I+diag(n,1)+diag(n,-1);
b = zeros(100,1);
x0 = zeros(100,1);
b(1) = 1;
b(100) = 1;
x0(1)=1;
TOL=1e-10;
timer1=tic;
error = Gauss_Seidel_iteration(A,b,x0);
fprintf('The consuming time of Gauss Seidel method = %f sec\n', toc(timer1)); 
timer2=tic;
error1 = jacobi_iteration(A,b,x0);
fprintf('The consuming time of Jacobi method = %f sec\n', toc(timer2));
x=[1:1e5];
semilogy(x,error,x,error1)
xlabel('iteration step')
ylabel('relative error')
legend('Gauss Seidel','Jacobi');
grid on
end