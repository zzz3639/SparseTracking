function [ tracksort, v ] = sorttrackid( trackori )
%Usage: [ tracksort, v ] = sorttrackid( trackori )
%   Sort a trajectory matrix [x,y,I,t,id] by id
%   for every trajectory(same id), locations are sorted by t
%   v is the order: tracksort = trackori(v,:)
[~,v1] = sort(trackori(:,4));
trackori = trackori(v1,:);
[~,v2] = sort(trackori(:,5));
tracksort = trackori(v2,:);
v = v1(v2,1);

end

