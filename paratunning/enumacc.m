function [ matchnum, matches ] = enumacc( trjtrue, trjansenum, thint )
%ENUMACC Summary of this function goes here
%   thint: dim particles under thint are removed

L = length(trjansenum);
matchnum = zeros(L,1);
matches = cell(L,2);

for i=1:L
    trjansthis = trjansenum{i};
    trjfilted = intensityfilter4column(trjansthis.track,thint);
    [m1,m2] = oneonematch(trjfilted,trjtrue.track);
    matchnum(i) = sum(m1>0);
    matches{i,1} = m1;
    matches{i,2} = m2;
end


end

