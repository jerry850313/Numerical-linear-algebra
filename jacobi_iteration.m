 function  [k,T] = jacobi_iteration(A,b,x0,tol,max1)
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
B = D\(L+U);
f = D\b;
x = B*x0+f;
T=norm(x-x0);
k = 1;
while T>=tol
    x0 = x;
    x = B*x0+f;
    T=norm(x-x0);
    k = k+1;
    if(k>=max1)
        disp(['It''s out of ',maxl,' times.'])
        return;
    end
end

x
end