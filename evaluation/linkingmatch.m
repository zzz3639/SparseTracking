function [m,s] = linkingmatch( trjtrue, trjans, th )
%LINKINGMATCH estimate the error rate of trjans based on trjtrue
%   trjtrue, trjans: trajectory matrix [x,y,I,t,id]
%   a link from emitter a to emitter b is a correct link if and only if a's
%     align in trjtrue is linked with b's align in trjtrue, and both
%     distance is not larger than th.

T = max(trjans(:,4));
trjtruesort = sorttrackid(trjtrue);
trjanssort = sorttrackid(trjans);
m = 0;
s = 0;
for t=0:T-1
    [idxanss,idxanst] = trjlinkcat(trjans,t);
    vt = find(trjtrue(:,4)==t);
    trjtrueframet = trjtrue(vt,:);
    vtp1 = find(trjtrue(:,4)==t+1);
    trjtrueframetp1 = trjtrue(vtp1,:);
    trjanslinkedt = trjans(idxanss,:);
    trjanslinkedtp1 = trjans(idxanst,:);
    [ alignresult_t ] = MolListAlign(trjanslinkedt(:,1:3), trjtrueframet(:,1:3));
    [ alignresult_tp1 ] = MolListAlign(trjanslinkedtp1(:,1:3), trjtrueframetp1(:,1:3));
    id1 = trjtrue(vt(alignresult_t(:,7)),5);
    D1 = alignresult_t(:,4);
    id2 = trjtrue(vtp1(alignresult_tp1(:,7)),5);
    D2 = alignresult_tp1(:,4);
    cl = (id1==id2&D1<th&D2<th);
    m = m+sum(cl);
    s = s+length(id1);
end


end

