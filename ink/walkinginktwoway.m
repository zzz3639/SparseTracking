function [ reswalking ] = walkinginktwoway( trj, knntrj, p0, q, m, tmax )
% random walk(monte carlo walk) for trajectory set trj start from position p0
% input:
%   trj: trajectories[x,y,I,t,id], no one frame trajectory allowed here.
%   knntrj: 
%   p0: start position
%   q: with probability q jump to k nearest neighbor, with probability 1-q
%      jump to the position of next (or previous, decided by current direction) frame
%   m: m monte carlo trials
%   tmax: monte carlo walking steps
% output:
%   reswalking: m*2*tmax;

% initialize
trj = sorttrackid( trj );
rng('shuffle');
kdtree = KDTreeSearcher(trj(:,1:2));
idx = knnsearch(kdtree,p0);
clear kdtree;
nmol = size(trj,1);
k = size(knntrj,2);

%distant matrix representing the molecules to their k nearest neighbor
Djump = zeros(nmol,k);
for i=1:k
    Djump(:,i) = sqrt(sum((trj(:,1:2) - trj(knntrj(:,i),1:2)).^2,2));
end

% in k nearest neighbor jump, the probability that jump to i'th nearest
% neighbor is proportional to 1/d_i
Pjump = repmat(1.0,nmol,k)./Djump;
Pjump = Pjump./repmat(sum(Pjump,2),1,k);

% add empty id for the judgement of wether or not a position is the start or end of a molecule.
trjidsearch = zeros(nmol+2,1);
trjidsearch(1) = trj(1,5)-1;
trjidsearch(end) = trj(end,5)+1;
trjidsearch(2:end-1,1) = trj(:,5);

p0start = trj(idx,1:2);
if trjidsearch(idx)~=trjidsearch(idx+1)
    directionthis = ones(m,1);
elseif trjidsearch(idx+1)~=trjidsearch(idx+2)
    directionthis = -ones(m,1);
else
    directionthis = 2*(randn(m,1)>0)-1;
end
directionnext = zeros(size(directionthis));

%monte carlo walkers in frame 1, m replicates
reswalking = zeros(m,2,tmax);
mcframethis = repmat(p0start,m,1);
reswalking(:,:,1) = mcframethis;
idxthis = repmat(idx,m,1);
idxnext = zeros(size(idxthis));
idxdirection = zeros(size(idxthis));
for i=2:tmax
    qsample = rand(m,1);
    %these walkers jump to nearest neighbor
    walkernn = find(qsample<q);
    %other walkers walk along this trajectory
    walkernormal = find(~(qsample<q));
    %normal walkers walk one step
    idxnext = idxthis;
    idxnext(walkernormal,1) = idxthis(walkernormal,1)+directionthis(walkernormal,1);
    %jump walkers jump to nearest neighbor
    jumpsample = mnrnd(1,Pjump(idxthis(walkernn,1),:));
    [jumpsituI,jumpsituJ,jumpsituS] = find(knntrj(idxthis(walkernn,1),:).*jumpsample);
    [u,v] = sort(jumpsituI);
    idxnext(walkernn,1) = jumpsituS(v);
    
    %define walking directions in the next step
    directionthis(walkernn,1) = 2*(randn(length(walkernn),1)>0)-1;
    directionnext = directionthis;
    
    %change direction if reaching the end
    idxdirection = idxnext + directionthis;
    v = find(trjidsearch(idxdirection+1,1)~=trjidsearch(idxnext+1,1));
    directionnext(v,1) = -directionnext(v,1);
    %record the positions of the walkers
    mcframethis = trj(idxnext,1:2);
    reswalking(:,:,i) = mcframethis;
    %update idxthis and directionthis
    idxthis = idxnext;
    directionthis = directionnext;
end

end



