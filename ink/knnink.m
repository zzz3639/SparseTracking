%{
function [ knntrj ] = knnink( trj, k )
%return k nearest neighbor index for all locations inside trj
%   inside this function, all non-self locations are valid nearest neighbor candidates

kdtree = KDTreeSearcher(trj(:,1:2));
idx = knnsearch(kdtree,trj(:,1:2),'K',k+1);
knntrj = idx(:,2:end);

end
%}

function [ knntrj ] = knnink( trj, k )
%return k nearest neighbor index for all locations inside trj
%   inside this function, all non-self locations are valid nearest neighbor candidates

N = size(trj,1);
kwidth = 2;
kdtree = KDTreeSearcher(trj(:,1:2));
searchset = [1:N]';
knntrj = zeros(N,k);
idtrjfull = trj(:,5);

while length(searchset)>0
    idx = knnsearch(kdtree,trj(searchset,1:2),'K',k*kwidth+1);
    idtrj = trj(searchset,5);
    idself = repmat(idtrj,1,k*kwidth);
    idnn = idtrjfull(idx(:,1:end-1));
    nonselfidx = firstknonzero( (idself~=idnn), k );

    Nsearch = length(searchset);
    knntrj(searchset,:) = catmatrix(idx, repmat([1:Nsearch]',1,k), nonselfidx);
    v = find(nonselfidx(:,end)==(k*kwidth+1));
    searchset = searchset(v,:);
    kwidth = kwidth+1;
end

end




