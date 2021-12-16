function [ tracksort, v ] = sorttrack( trackori )
%Usage: [ tracksort ] = sorttrack( trackori )
%   Sort a trajectory matrix [x,y,I,t,id] by t
%   v is the order: tracksort = trackori(v,:)
[~,v1] = sort(trackori(:,5));
trackori = trackori(v1,:);
[~,v2] = sort(trackori(:,4));
tracksort = trackori(v2,:);
v = v1(v2,1);

end

