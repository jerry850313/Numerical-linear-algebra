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
        t=t+1;
        x0=x;
        if t>(1e5-1)
            return
        end
    end
end