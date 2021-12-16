function [ trjlink ] = linktwotrj( trj1, trj2 )
% trj1 and trj2 are two trajectory lists, generate trjlink, where frame
% 0-t1 same as trj1, t1+1 to t2 same as trj2
%   trj1 = [x,y,I,t,id], 0<=t(i)<=t1
%   trj2 = [x,y,I,t,id,link], t1+1<=t(i)<=t2
%      i'th(line number) emitter in trj2 link to link(i)'th
%      emitter in trj1, link(i)=0 means i'th molecule in trj2 links to no one
%   trjlink = [x,y,I,t,id], trajectories after linking

% find non-zero values of link
v = find(trj2(:,6)>0);
% verify link (one id, one link)
linkunique = unique(trj2(v,6));
if length(v)~=length(linkunique)
    trjlink = zeros(0,5);
    return;
end

% assign id to non-zero link elements (trj1 link to trj2)
idmap = containers.Map('KeyType', 'int64', 'ValueType', 'int64');
for i=1:length(v)
    idmap(trj2(v(i),5)) = trj1(trj2(v(i),6),5);
end

% find all ids for trj2
idtrj2 = unique(trj2(:,5));
% assign id to zero link elements (trajectories in trj2 with no predecessor)
idmax = max(trj1(:,5));
k=1;
idtrj2num = length(idtrj2);
for i=1:idtrj2num
    idthis = idtrj2(i);
    if idmap.isKey(idthis)
    else
        idmap(idthis) = idmax + k;
        k = k + 1;
    end
end

% write new ids to trj2 
trj2linenum = size(trj2,1);
for i=1:trj2linenum
    trj2(i,5) = idmap(trj2(i,5));
end

% union trj1 and trj2
trjlink = [trj1;trj2(:,1:5)];

end

