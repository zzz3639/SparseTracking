function [ mcmv ] = resconfinedprocessing( trjans, S, startpoint, Tmaxmc, stride, pjump, knnjump, model, corenum )
% Usage: [ mcmv ] = resconfinedprocessing( trjans, startpoint, Tmaxmc, stride, pjump, corenum )
%  Run inking algorithm based on trjans from startpoint.
%  Inputs:
%   trjans: trajectories [x,y,I,t,id]
%   S: size of movie board, [S(1),S(2)]
%   startpoint: where the monte carlo particles start from
%   Tmaxmc: how many steps monte carlo would run
%   stride: We don't have to record every monte carlo step, record
%      (Tmaxmc/stride) steps
%   pjump: probability monte carlo particle jump to knn in one step
%   knnjump: k nearest neighbor k
%   model: 'one' for walkinginkoneway, 'two' for walkinginktwoway
%   corenum: number of computer cores to run this, large number takes a
%      lot of memory
%  Outputs:
%   mcmv: the resulting movie, size: S(1)*zoom,S(2)*zoom,(Tmaxmc/stride);
zoom = 10;
mparticle = 20000;
Lfilter = 3;

[ reswalking ] = runwalkingink( trjans, knnjump, Lfilter, startpoint, pjump, mparticle, Tmaxmc, model, corenum );
[ trjmc ] = reswalking2trj( reswalking );
clear reswalking;
[ mcmv ] = visualizetrjtrack( trjmc, S, 0, zoom ,stride);
clear trjmc;
m = 4;
for i=1:m
    fprintf('\nouterloop: i=%d\n',i);
    [ reswalking ] = runwalkingink( trjans, knnjump, Lfilter, startpoint, pjump, mparticle, Tmaxmc, model, corenum );
    [ trjmc ] = reswalking2trj( reswalking );
    clear reswalking;
    [ mcmvtemp ] = visualizetrjtrack( trjmc, S, 0, zoom ,stride);
    clear trjmc;
    mcmv = mcmv + mcmvtemp;
end


end

