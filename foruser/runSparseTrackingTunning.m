function [ trjans, trjansdetail ] = runSparseTrackingTunning( moviethis, movietunning, trjtruetunning, sig, psfdecay, bsize, D, alpha, thint, corenum, jobname )
% Run our sparse tracking algorithm, hyper-parameters are tunned inside
% input:
%   moviethis: image series to be processed
%   movietunning: image series for hyper-parameter tunning
%   trjtruetunning: ground truth for movietunning
%   sig: psf width, gaussian sigma
%   psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%   bsize: boundary size
%   D: diffusion rate Gaussian sigma
%   alpha: intensity linkage
%   thint: dim particles under thint are removed
%   corenum: number of computer cores for parfor
%   jobname: a string, as tracksolverlspfull writes and reads files from hard disk, 
%            we use jobname to ensure 2 jobs use different files, avoid conflict.
% outputs:
%   trjans: molecule locations in all the frames,
%        trjans.no is T*1 array records noise values of all frames
%        trjans.track is a matrix, [x,y,I,T], 2 position values, intensities, frame
%   trjansdetail: record the detailed information of sequential on off loop
%

[ spratebest, lambdaAbest, lambdaBbest, matchnum, trjansenum, trjans5cenum ] = runparatunning( movietunning, trjtruetunning, sig, psfdecay, bsize, D, alpha, thint, corenum );
[ trjans, trjansdetail ] = runtrackingfixpara( moviethis, sig, psfdecay, bsize, D, alpha, spratebest, lambdaAbest, lambdaBbest, corenum, jobname );

[ ~, ~, ~, ~, ~, ~, bcmax ] = defaultparaset();
trjans = trjansdetailtotrjans(trjansdetail, bcmax, thint);

end

