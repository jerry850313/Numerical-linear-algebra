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