function [ trjremain ] = rmshorttrj( trjfull, l )
% trajectories shorter than l (trjlen<l) are removed

[trjunique,lbnum]= molidunique(trjfull);
trjsort = sorttrackid(trjunique);
L = zeros(lbnum,1);
rmid = zeros(size(trjsort,1),1);
for i=0:lbnum-1
    v = find(trjsort(:,5)==i);
    L(i+1) = length(v);
    if L(i+1)<l
        rmid(v,1) = 1;
    end
end 

v = find(rmid==0);

trjremain = trjsort(v,:);
[trjunique,lbnum]= molidunique(trjremain);
trjremain = sorttrackid(trjunique);

end

