function [m,s] = linkingmatchkframe( trjtrue, trjans, thdist, K )
% LINKINGMATCH estimate k-frame linking error rate of trjans based on trjtrue
% k-frame link consists of k positions, k-1 links
%   trjtrue, trjans: trajectory matrix [x,y,I,t,id]
%   a link from emitter a to emitter b is a correct link if and only if a's
%     align in trjtrue is linked with b's align in trjtrue, and both
%     distance is not larger than th.

T = max(trjans(:,4));
% sort by id
[trjtruesort, ordertrue] = sorttrackid(trjtrue);
[trjanssort, orderans] = sorttrackid(trjans);

[ linkinfo, orderlink ] = trjlinkinfo( trjanssort );
% reorder
trjanssort = trjanssort(orderlink,:);
% align to ground truth
[ linealigned ] = trjaligntotrue( trjtruesort, trjanssort );
% compute alignment distance
D = Dist2D(trjanssort(:,1:2), trjtruesort(linealigned,1:2));
idfromtrue = trjtruesort(linealigned,5);

% find valid align which are closer than thdist
alignvalid = (D<thdist);

% find all k-frame links
[linkkframe] = findlinkkframe(linkinfo(:,1),K);

% compute linking error
% distance validation
[alignkvalid] = findalignkvalid(alignvalid,K);
% id validation
[idkvalid] = findidkvalid(idfromtrue,K);
% number of correct k-frame links
m = sum(linkkframe.*alignkvalid.*idkvalid);
% number of k-frame links
s = sum(linkkframe);

end

function [D] = Dist2D(P1,P2)
D = sum((P1-P2).*(P1-P2),2);
D = sqrt(D);
end

function [linkk] = findlinkkframe(linkstates,K)
% find k-frame links based on links
K = K-1;
if size(linkstates,1)==1
    linkstates = linkstates';
end
N = length(linkstates);
if N<K
    linkk = zeros(N,1);
    return;
end
linksum = linkstates;
linktemp = linkstates;
for i=2:K
    linktempp1 = [linktemp;0];
    linktemp = linktempp1(2:end,1);
    linksum = linksum + linktemp;
end

linkk = (sum(linksum,2)==K);

end

function [alignkvalid] = findalignkvalid(alignvalid,K)
% Compute validation of k-fames links based on aligned distance, mathematically 
% equal to the problem of finding k-frame links
[alignkvalid] = findlinkkframe(alignvalid,K+1);
end

function [idkvalid] = findidkvalid(idfromtrue,K)
% Compute validation of k-fames links based on id: k positions must align to the same emitter in trjtrue
% mathematically equal to the problem of finding k-frame links
if size(idfromtrue,1)==1
    idfromtrue = idfromtrue';
end
idp1 = [idfromtrue(2:end,1);idfromtrue(end)+1];
idsame = (idfromtrue==idp1);
[idkvalid] = findlinkkframe(idsame,K);

end



