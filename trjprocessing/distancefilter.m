function [ trjfilted, lbnum ] = distancefilter( trj, thdist )
% links longer than thdist are removed
%   inputs
%     trj: [x,y,I,t,id], trajectories to be filterd
%     thint: non-negative number
%   outputs:
%     trjfilted: [x,y,I,t,id], trajectories after the filter
%     lbnum: number of ids in trjfilted

% sort trj by id
[trj] = moliduniquesorted(trj);
% get linking information
[ linkinfo, order ] = trjlinkinfo( trj );
trj = trj(order,:);
linkstates = linkinfo(:,1);
D = linkinfo(:,4);
% find long linkings
cut = ((linkstates==1)&(D>thdist));
cut = [0;cut(1:end-1,1)];
% remove dim localizations and reassign ids
newid = cumsum(cut) + trj(:,5);
trj(:,5) = newid;
[trjfilted,lbnum] = moliduniquesorted(trj);

end
