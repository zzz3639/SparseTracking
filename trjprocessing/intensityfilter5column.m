function [ trjfilted, lbnum ] = intensityfilter5column( trj, thint )
% emitters whose intensity lower than thint are removed, trajectories with 
%   particles removed in the middle are cut into segments.
%   inputs
%     trj: [x,y,I,t,id], trajectories to be filterd
%     thint: non-negative number
%   outputs:
%     trjfilted: [x,y,I,t,id], trajectories after the filter
%     lbnum: number of ids in trjfilted

% sort trj by id
[trj] = moliduniquesorted(trj);
% find dim localizations
u = (trj(:,3)<thint); % dim lines
linenotdim = find(~u);
% remove dim localizations and reassign ids
newid = cumsum(u) + trj(:,5);
trj(:,5) = newid;
trjfilted = trj(linenotdim,:);
[trjfilted,lbnum] = moliduniquesorted(trjfilted);

end
