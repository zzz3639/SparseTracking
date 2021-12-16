function [ nndist ] = nndistribution( trj )
% Compute nearest neighbor distribution
%  Usage:
%    [nndist] = nndistribution( trj )
%    trj: trajectory matrix, [x,y,I,t,id]
%    nndist: nearest neighbor distribution

trjsort = sorttrack( trj );
nnmax = 100;

framelist = unique(trjsort(:,4));
N = length(framelist);
nnvalue = zeros(size(trj,1),1);
k=1;
for i=1:N
    t = framelist(i);
    v = find(trjsort(:,4)==t);
    l = length(v);
    if l==1
        nnvalue(k) = nnmax;
        k = k+1;
        continue;
    end
    trjframe = trjsort(v,:);
    D = dist(trjframe(:,1:2),trjframe(:,1:2));
    D(logical(eye(size(D)))) = 100;
    nnthis = min(D);
    nnvalue(k:k+l-1,1) = nnthis';
    k = k+l;
end

nndist = nnvalue;

end


function D=dist(mu1,mu2)
    m=size(mu1,1);
    n=size(mu2,1);
    D=zeros(n,m);
    for i=1:m
        T=mu2-repmat(mu1(i,:),n,1);
        D(:,i)=T(:,1).*T(:,1)+T(:,2).*T(:,2);
    end
    D=sqrt(D);
end

