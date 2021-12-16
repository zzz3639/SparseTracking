function [idxs,idxt] = trjlinkcat(trj,ts)
% find all emitters in frame ts who have link to ts+1, emitters in ts is
% recorded in idxs, emitters in ts+1 is recorded in idxt

vs = find(trj(:,4)==ts);
vt = find(trj(:,4)==(ts+1));
ids = trj(vs,5);
idt = trj(vt,5);
vc = [vs;vt];
idcforsort = [ids;idt+0.1];
[sortedid,vid] = sort(idcforsort);
pairid = sortedid(2:end,:) - sortedid(1:end-1,:);
pairedline = find(pairid<0.2);
idxs = vc(vid(pairedline,:));
idxt = vc(vid(pairedline+1,:));

end