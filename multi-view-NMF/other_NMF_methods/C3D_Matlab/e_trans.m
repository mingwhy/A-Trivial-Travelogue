function [V,Lamda]=e_trans(D)
% Compute S and the corresponding eigensystems for the symmetric(co-expression) input data matrices D
% Input:
%     D:    each cell D{i} is a symmetric co-expression matrix
% Output:
%     V:    the common factor of hosgvd decomposition of D, 
%     each column is an eigenvector   
%     Lamda: each diagonal value is an eignevalue 
% W is computed directly from the D
%
% Written by Xiaolin Xiao, 2012


if nargin < 1 || isempty(D)
    error(message('TooFewInputs'));
end

cond=length(D);
if cond<2
    error('number of conditions is less than 2');
end

drsize=nan(1,cond);
real_tag=nan(1,cond);
for i=1:cond
    drsize(i)=size(D{i},1);
    real_tag(i)=isreal(D{i});
end
if isequal(drsize, drsize(1).*ones(1,cond))==0
    error(message('MATLAB:hogsvd:MatrixRowMismatch'))
end
if isequal(real_tag,ones(1,cond))==0
    error('Complex numbers exist');
end


A=D; % take co-expression

S=zeros(length(A{1}));
for i=1:cond-1
    for j=i+1:cond
        if size(D{j},1)<size(D{j},2)
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*fast_pinv(D{j}')*fast_pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        else
            if size(D{i},1)<size(D{i},2)
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*fast_pinv(D{i}')*fast_pinv(D{i}));
            else
                S=S+(A{i}*pinv(D{j}')*pinv(D{j})+A{j}*pinv(D{i}')*pinv(D{i}));
            end
        end           
    end
end
S=S./(cond*(cond-1));
rs=rank(S);

if rs<length(S)-1
    [V,Lamda]=eigs(S,rs);
else
    [V,Lamda]=eig(S);
end

if isreal(Lamda)==0
    lam=diag(Lamda);
    k=1;
    while isreal(lam(k))==1
        k=k+1;
    end
    if k<=1
        error('non real patterns found');
    else
        Lamda=Lamda(1:k-1,1:k-1);
        V=V(:,1:k-1);
    end
end


colnorm=sqrt(sum(V.^2,1)); 
colnorm=repmat(colnorm,[size(V,1) 1]);
V=V./colnorm; 


end % end of function








