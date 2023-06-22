function Y = fast_pinv(G)
% Fast computation of moore-penrose inverse of matrix G(m*n) if m<n
% reference: P. Courrieu, Fast Computation of Moore-Penrose Inverse
% Matrices, Neural Information Processing -Letters and Reviews, 8 (2),
% 2005.


[m,n]=size(G); 
transpose=false;
if m<n
    transpose=true;
    A=G*G';
    n=m;
else
    A=G'*G;
end

dA=diag(A); 
tol= min(dA(dA>0))*1e-9;
L=zeros(size(A));
r=0;
for k=1:n
    r=r+1;
    L(k:n,r)=A(k:n,k)-L(k:n,1:(r-1))*L(k,1:(r-1))';
    if L(k,r)>tol
        L(k,r)=sqrt(L(k,r));
        if k<n
            L((k+1):n,r)=L((k+1):n,r)/L(k,r);
        end
    else
        r=r-1;
    end
end
L=L(:,1:r);

M=inv(L'*L);
if transpose
    Y=G'*L*M*M*L';
else
    Y=L*M*M*L'*G';
end