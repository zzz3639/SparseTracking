function [trackidunique,lbnum]= molidunique(trackori)
% Usage: [trackidunique,lbnum]= molidunique(trackori)
% this function unique the molecule ids and re-assign them to 0:lbnum-1
% inputs:
%  trackori: [x,y,I,t,id]
% outputs:
%  trackidunique: [x,y,I,t,id]
%  lbnum, number of ids
    trackidunique = trackori;
    lbindex = unique(trackori(:,5));
    lbnum = length(lbindex);
    k = 0;
    for i=1:lbnum
        u = find(trackori(:,5)==lbindex(i));
        trackidunique(u,5) = i-1;
    end
    trackidunique = sorttrack(trackidunique);
end