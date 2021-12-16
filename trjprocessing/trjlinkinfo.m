function [ linkinfo, order ] = trjlinkinfo( trj )
% return the linking information of a trajectory matrix
%   input:
%     trj: trajectory matrix [x,y,I,t,id]
%   output:
%     linkinfo: [link,dx,dy,d], line number same to trjsorted, link(i)=0/1,
%       recording whether trjsorted line i link to i+1
%     order: trj was sorted by id inside, this is the sort order.

[trjsorted, order] = sorttrackid(trj);
N = size(trjsorted,1);
linkinfo = zeros(N,4);

Pnext = [trjsorted(2:end,1:2);zeros(1,2)];
D = Dist2D(trjsorted(:,1:2),Pnext);
linkinfo(:,4) = D;
linkinfo(:,2:3) = Pnext-trjsorted(:,1:2);

id = trjsorted(:,5);
idnext = [id(2:end,1);id(end,1)+1];
linkinfo(:,1) = (id==idnext);

end

function [D] = Dist2D(P1,P2)
D = sum((P1-P2).*(P1-P2),2);
D = sqrt(D);
end
