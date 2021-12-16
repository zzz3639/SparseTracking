function [ trjans ] = trjansdetailtotrjans( trjansdetail, bcmax, thfilter )
%for a movie with n frames, trjansdetail is a cell array recording the
%tracking results of submovie (1->bcmax+1), (2->bcmax+2), (3->bcmax+3) ...
%  inputs:
%   trjansdetail: cell array, trjansdetail{i}.track = [x,y,I,t,id],trjansdetail{i}.no=[no1;no2,...]
%   bcmax: bcmax+1 is the length of submovies
%   thfilter: non-negative number, emitters whose intensity lower than thint are removed
%  outputs:
%   trjans: trjans.track=[x,y,I,t,id], trjans.no=[no1;no2...] 

% default parameters
alphaalign = log(2); % to align trjansdetail{i} and trjansdetail{i+1}
Dalignrate = 1/3;
thalign = 1.0;

N = length(trjansdetail);
if N==1
    trjans = trjansdetail{1};
    return;
end
% get the information
[Dave,Ifave] = trjansdetailinfo(trjansdetail, thfilter);
Dalign = Dalignrate*Dave;

% sort and filter every element in trjansdetail by id
for i=1:N
    trjansdetail{i}.track = sorttrackid(trjansdetail{i}.track);
    trjansdetail{i}.track = intensityfilter5column( trjansdetail{i}.track, thfilter );
end
% 1--mdinter are defined by the first element
Tmax = bcmax+N;
trjans.track = zeros(0,5);
trjans.no = zeros(0,1);
mdinter = floor((2+bcmax)/2);
trjtemp = trjansdetail{1};
v = find(trjtemp.track(:,4)<mdinter);
trjans.track = [trjans.track;trjtemp.track(v,:)];
trjans.no = [trjans.no;trjtemp.no(1:mdinter,1)];
t = mdinter-1; %we have linked 0-->mdinter-1 
% relink inside this loop
for i=2:N
    % emitters in frame t in trjans
    vto = find(trjans.track(:,4)==t);
    trjtemp = trjansdetail{i};
    % emitters in trjansdetail{i} from frame mdinter-2 to mdinter-1
    % vfromnext correspond to the next frame (t+1).
    [vfrom,vfromnext] = trjlinkcat(trjtemp.track,mdinter-2);
    %emitter positions and intensities
    emitterto = trjans.track(vto,1:3);
    emitterfrom = trjtemp.track(vfrom,1:3);
    %frame t in trjans correspond to frame mdinter-2 in trjtemp
    [alignfrom,alignto] = emitteralias(emitterfrom,emitterto,Dalign,alphaalign,thalign);
    [alignfrom,alignto] = checkmapping(alignfrom,alignto);
    vfromnextfull = find(trjtemp.track(:,4)==mdinter-1);
    vfromlink = subsetsearch(vfromnextfull,vfromnext);
    trjnextframe = trjtemp.track(vfromnextfull,:);
    trjnextframe(:,4) = t+1;
    link = zeros(length(vfromnextfull),1);
    link(vfromlink(alignfrom)) = vto(alignto);
    trjnextframe = [trjnextframe,link];
    trjans.track = linktwotrj( trjans.track, trjnextframe );
    trjans.no = [trjans.no;trjtemp.no(mdinter,1)];
    t = t+1;
end

% link the remaining frames
    % emitters in frame t in trjans
vto = find(trjans.track(:,4)==t);
    % emitters in trjansdetail{end} from frame mdinter-1 to mdinter
    % vfromnext correspond to the next frame (t+1).
trjtemp = trjansdetail{end};
[vfrom,vfromnext] = trjlinkcat(trjtemp.track,mdinter-1);
    %emitter positions and intensities
emitterto = trjans.track(vto,1:3);
emitterfrom = trjtemp.track(vfrom,1:3);
    %frame t in trjans correspond to frame mdinter-2 in trjtemp
[alignfrom,alignto] = emitteralias(emitterfrom,emitterto,Dalign,alphaalign,thalign);
[alignfrom,alignto] = checkmapping(alignfrom,alignto);
    %the remaining frames
vremain = find(trjtemp.track(:,4)>=mdinter);
vfromlink = subsetsearch(vremain,vfromnext);
trjnextframe = trjtemp.track(vremain,:);
trjnextframe(:,4) = trjnextframe(:,4)-bcmax+Tmax-1;
    %the link
link = zeros(length(vremain),1);
link(vfromlink(alignfrom)) = vto(alignto);
trjnextframe = [trjnextframe,link];
trjans.track = linktwotrj( trjans.track, trjnextframe );
trjans.no = [trjans.no;trjtemp.no(mdinter+1:end,1)];

end


function [idline_s,idline_t] = idlookup(ids)
%ids is n*1 sorted vector
%idline_s: each id start from which line
%idline_t: each id end at which line
%example: idlookup([1;1;2;2;3;3]) returns [1,3,5] and [2,4,6], 
%  label 1: ids(1:2)
%  label 2: ids(3:4)
%  label 3: ids(5:6)
    ids1 = [ids(1)-1;ids];
    ids2 = [ids;ids(end)+1];
    v = find(ids1~=ids2);
    idline_s = v(1:end-1,1);
    idline_t = v(2:end,1)-1;
end
