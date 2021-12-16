function [D,If,SVarray] = trjinfo( trj, noprint )
%show the summary of a trajectory matrix
%   trj: [x,y,I,t,id]

trj = sorttrackid(trj);
[trj,lbnum]= moliduniquesorted(trj);
%number of frames
Nframe = max(trj(:,4))+1;
if ~exist('noprint')
    fprintf('\n Number of frames: %d',Nframe);
end

%number of molecules
Nmol = lbnum;
if ~exist('noprint')
    fprintf('\n Number of molecules: %d',Nmol);
end

%number of links
Nlk = size(trj,1) - lbnum;
if ~exist('noprint')
    fprintf('\n Number of links: %d',Nlk);
end

%average survival length of the molecules
[idline_s,idline_t] = idlookup(trj(:,5));
SVarray = (idline_t-idline_s)+1;
SVlen = size(trj,1)/Nmol;
if ~exist('noprint')
    fprintf('\n Average survival length: %f',SVlen);
end
%average intensity of the molecules
Intave = mean(trj(:,3));
if ~exist('noprint')
    fprintf('\n Average intensity of the molecules: %f',Intave);
end

% separate the trajectories of the molecules
lbtrjfrom = [trj(1,5)-1;trj(1:end-1,5)];
lbtrjto = trj(:,5);
lbsame = (lbtrjfrom==lbtrjto);
vsame = find(lbsame);

% Displacement
Ltrjfrom = [trj(1,1:2);trj(1:end-1,1:2)];
Ltrjto = trj(:,1:2);

Dtrj = dist(Ltrjfrom,Ltrjto);
Dave = sqrt(sum(Dtrj.*Dtrj.*lbsame)/sum(lbsame));
D = Dtrj(vsame,:);
if ~exist('noprint')
    fprintf('\n Average squared displacement: %f', Dave);
    fprintf('\n Median displacement: %f', median(D));
end

% Intensity fluctuation
Iftrjfrom = [trj(1,3);trj(1:end-1,3)];
Iftrjto = trj(:,3);
Ifdist = abs(log(Iftrjto) - log(Iftrjfrom));
If = Ifdist(vsame,:);
Ifave = sum(Ifdist.*lbsame)/sum(lbsame);
if ~exist('noprint')
    fprintf('\n Intensity fluctuation: %f', Ifave);
end


if ~exist('noprint')
    fprintf('\n');
end

end

function D=dist(L1,L2)
%Euclid distance, L1,L2: n*2 
    D = (L2(:,1)-L1(:,1)).*(L2(:,1)-L1(:,1)) + (L2(:,2)-L1(:,2)).*(L2(:,2)-L1(:,2));
    D = sqrt(D);
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


