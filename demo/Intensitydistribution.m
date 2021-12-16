function [] = Intensitydistribution( trj )
%INTENSITYDISTRIBUTION Summary of this function goes here
%   trj: [x,y,I,t,id]

[trjunique,lbnum]= molidunique(trj);
trjsort = sorttrackid(trjunique);
N = size(trjsort,1);
Int = trjsort(:,3);
m = mean(log(Int));
D = var(log(Int));


L = zeros(lbnum,1);
for i=0:lbnum-1
    L(i+1) = length(find(trjsort(:,5)==i));
end
Lf = N-lbnum;
I1 = zeros(Lf,1);
I2 = zeros(Lf,1);
k=1;
l=1;
for i=1:lbnum
    for j=1:L(i)-1
        I1(l) = trjsort(k,3);
        I2(l) = trjsort(k+1,3);
        k=k+1;
        l=l+1;
    end
    k=k+1;
end
[us,vs] = sort(L);
clear us;
numlentrack=20;
Ilengthmax = cell(numlentrack,1);
for i=1:numlentrack
    idlengthmax = vs(end-i+1)-1;
    v = find(trjsort(:,5)==idlengthmax);
    Ilengthmax{i} = trjsort(v,3);
end


end



