function HW_2_2
clc;clear all;
TOL=1e-10;
M=3e5;
function  [k,T] = Gauss_Seidel_iteration(A,b,x0,tol,max1)
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
Tg = -(D+L)\U;
Cg = (D+L)\b;
x1=x0;
x = Tg*x0+Cg;
T=norm(x-x0);
t=norm(x0);
k = 1;
while (T/t)>=tol
    x0 = x;
    x = Tg*x0+Cg;
    T=norm(x-x0);
    t=norm(x0);
    k = k+1;
    if(k>=max1)
        return
    end
    T=norm(x-x1);
end
end
 function  [k,T] = jacobi_iteration(A,b,x0,tol,max1)
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
B = D\(L+U);
f = D\b;
x1=x0;
x = B*x0+f;
T=norm(x-x0);
t=norm(x0);
k = 1;
while (T/t)>=tol
    x0 = x;
    x = B*x0+f;
    T=norm(x-x0);
    t=norm(x0);
    k = k+1;
    if(k>=max1)
        return
    end
    T=norm(x-x1);
end
end
v=1:0.2:2;
A = fliplr(vander(v));
x0=zeros(6,1);
x0(1)=1;
b=A(1:6,1)+A(1:6,2)+A(1:6,3)+A(1:6,4)+A(1:6,5)+A(1:6,6);
[k1,error1]=jacobi_iteration(A,b,x0,TOL,M);
fprintf('Jacobi method uses %d times and its error is %d \n',k1,error1);
[k2,error2]=Gauss_Seidel_iteration(A,b,x0,TOL,M);
fprintf('Gauss Seidel method uses %d times and its error is %d \n',k2,error2);
end
