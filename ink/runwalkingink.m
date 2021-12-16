function [ reswalking ] = runwalkingink( trj, k, l, p0, q, m, tmax, model, corenum )
%RUNWALKINGINK Summary of this function goes here
%   k: k-nearest neighbor   
%   l: trajectories shorter than l (trjlen<l) are removed
%   model: 'one' for walkinginkoneway, 'two' for walkinginktwoway, 
%     'dance' for walkinginkdancing

trjremain = rmshorttrj( trj, l );
clear trj;
knntrj = knnink( trjremain, k );
mpara = zeros(corenum,1);
mpara = parajobnum( m, corenum );
jobt = cumsum(mpara);
jobs = [0;jobt(1:end-1,1)]+1;
reswalkingcell = cell(corenum,1);

%48 hours for parallel workers.
TimeOut=2880;
parpool('local',corenum,'IdleTimeout',TimeOut);

parfor i=1:corenum
    i
    if strcmp(model,'one')
        reswalkingcell{i} = walkinginkoneway( trjremain, knntrj, p0, q, mpara(i), tmax );
    elseif strcmp(model,'two')
        reswalkingcell{i} = walkinginktwoway( trjremain, knntrj, p0, q, mpara(i), tmax );
    else
        reswalkingcell{i} = walkinginkdancing( trjremain, knntrj, p0, q, mpara(i), tmax );
    end
end

%close the workers
delete(gcp('nocreate'));

reswalking = zeros(m,2,tmax);
for i=1:corenum
    reswalking(jobs(i):jobt(i),:,:) = reswalkingcell{i};
end
clear reswalkingcell;

end

