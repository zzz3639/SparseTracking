function [ linkingmatch ] = linkingmatchtrue( trj1, trj2)
% estimated linking error for falcon results(results with localizations only)
%   trj1: [x,y,I,t,id]; true trajectory answer
%   trj2 should have at least 4 columns [x,y,I,t];
%   linkingmatch: array of length size(trj1,1)-1. 
%     linkingmatch(i) recording wether or not trj1(i,:) to trj1(i+1,:) is a correct linking
%       yes if linkingmatch(i)==1
%       no if linkingmatch(i)==0
%       not even belong to the same trajectory if linkingmatch(i)==-1

trj1 = sorttrackid(trj1);
[m1,m2] = oneonematch(trj1,trj2);
n = length(m1);
linkingmatch = zeros(n-1,1);
for i=1:n-1
    if trj1(i+1,5)~=trj1(i,5)
        linkingmatch(i) = -1;
    elseif m1(i+1)>0&&m1(i)>0
        linkingmatch(i) = 1;
    else
        linkingmatch(i) = 0;
    end
end


end

