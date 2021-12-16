function [ spbest ] = falconparaselect( movietunning, trjtruetunning, thint )
%FALCONPARASELECT Summary of this function goes here
%   Detailed explanation goes here
sparray = falconpara();
N = length(sparray);
falcontunningcell = cell(N,1);
matches = zeros(N,3);
matches(:,2) = repmat(size(trjtruetunning.track,1),N,1);
for i=1:N
    sp = sparray(i);
    [molresults] = runtrackingfalcon(movietunning,sp);
    matches(i,3) = size(molresults,1);
    molresfilted = intensityfilter4column(molresults,thint);
    falcontunningcell{i} = molresfilted;
    [m1,m2] = oneonematch(trjtruetunning.track,molresfilted);
    matches(i,1) = sum(m1>0);
end

Acc = matches(:,1)./(matches(:,2)+matches(:,3));
[~,v] = max(Acc);
spbest = sparray(v);

end

