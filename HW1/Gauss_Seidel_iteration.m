function  [k T] = Gauss_Seidel_iteration(A,b,x0,tol,max1)
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
Tg = -(D+L)\U;
Cg = (D+L)\b;
x = Tg*x0+Cg;
T=norm(x-x0);
t=norm(x);
k = 1;
while (T/t)>=tol
    x0 = x;
    x = Tg*x0+Cg;
    T=norm(x-x0);
    t=norm(x);
    k = k+1;
    if(k>=max1)
        disp(['It''s out of ',maxl,' times.'])
        return;
    end
end
x
end