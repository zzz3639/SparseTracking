function [ trjans, trjansdetail ] = runtrackingfixpara( moviethis, sig, psfdecay, bsize, D, alpha, sprate, lambdaA, lambdaB, corenum, jobname )
% Using fixed parameters to run tracking algorithm
%  inputs:
%    moviethis: movie to be solved, m*n*T
%    sig: Gaussian sigma
%    psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%    bsize: boundary size
%    D: diffusion rate Gaussian sigma
%    alpha: intensity linkage
%    corenum: number of cores using in this program.
%    jobname: a string, as tracksolverlspfull writes and reads files from hard disk, 
%             we use jobname to ensure 2 jobs use different files, avoid conflict.

thfilter = 1; % emitter intensity lower than this are removed
Tmax = size(moviethis,3);
[ spratetunning, lambdaAtunning, lambdaBtunning, deltaIvalue, itemax, Kextend, bcmax ] = defaultparaset();
corenumideal = max(1,Tmax-bcmax);
corenum = min(corenum,corenumideal);
clear spratetunning;
clear lambdaAtunning;
clear lambdaBtunning;
trjtrue = emptytrack();
% prepare for parallele assignment
jobnum = Tmax-bcmax;
u = floor(jobnum/corenum);
umod = mod(jobnum,corenum);
barrel = u*ones(corenum,1); % how many jobs each core take
barrel(1:umod,1) = barrel(1:umod,1)+1;
jobcut = cumsum(barrel);
jobs = [0;jobcut(1:end-1,:)]+1;
jobt = jobcut;
trjansdetailcell = cell(corenum,1);
lambda = [lambdaA + [0:Tmax-1] * lambdaA*lambdaB, 0];
deltaI = deltaIvalue*ones(1,Tmax);
%72 hours for parallel workers.
TimeOut = 4320;
parpool('local',corenum,'IdleTimeout',TimeOut);
%run parallel computation
parfor i=1:corenum
    s = jobs(i);
    t = jobt(i)+bcmax;
    t = min(t,Tmax);
    movietask = moviethis(:,:,s:t);
    [ ~, ~, trjansdetailtask ] = sequentialonoff( movietask, sig, psfdecay, bsize, D, alpha, lambda, deltaI, sprate, itemax, Kextend, bcmax, trjtrue, 'no', [jobname,num2str(i)] );
    trjansdetailcell{i} = trjansdetailtask;
end

delete(gcp('nocreate'));

% concatenate the trjansdetail cell array
trjansdetail = cell(0,1);
for i=1:corenum
    trjansdetail = [trjansdetail;trjansdetailcell{i}];
end

%trjans = 'Use function trjansdetailtotrjans to reconstruct trjans';
trjans = trjansdetailtotrjans(trjansdetail, bcmax, thfilter);

end

