function [trackidunique,lbnum]= moliduniquesorted(trackori)
% Usage: [trackidunique,lbnum]= molidunique(trackori)
% this function unique the molecule ids and re-assign them to 0:lbnum-1
%  similar to function molidunique, faster, trackori will be sorted by id
% inputs:
%  trackori: [x,y,I,t,id]
% outputs:
%  trackidunique: [x,y,I,t,id]
%  lbnum, number of ids

trackidunique = sorttrackid(trackori);
[idline_s,idline_t] = idlookup(trackidunique(:,5));
lbnum = length(idline_s);
for i=1:lbnum
    trackidunique(idline_s(i):idline_t(i),5) = i-1;
end

end

function [idline_s,idline_t] = idlookup(ids)
%ids is n*1 sorted vector
    ids1 = [ids(1)-1;ids];
    ids2 = [ids;ids(end)+1];
    v = find(ids1~=ids2);
    idline_s = v(1:end-1,1);
    idline_t = v(2:end,1)-1;
end