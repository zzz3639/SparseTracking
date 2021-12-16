function [ srange,trange ] = molranget( track )
%Given trajectory matrix, molranget reports the start and end frame of each molecule
% input:
%   track: trajectory matrix, [x,y,I,t,id], t and id start at 0
% output:
%   srange: where the molecules start, idnum*1 array
%   trange: where the molecules end, idnum*1 array

maxtrack = sparse(track(:,4)+1,track(:,5)+1,track(:,4));
exitrack = sparse(track(:,4)+1,track(:,5)+1,ones(size(track,1),1));
Lrange = [full(sum(exitrack,1))]';
trange = [full(max(maxtrack,[],1))]';
srange = [trange-Lrange+1];

end

