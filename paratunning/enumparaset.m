function [ trjansenum, mvfullenum, trjansdetailenum ] = enumparaset( moviethis, trjtrue, sig, psfdecay, bsize, D, alpha, sprate, lambdaA, lambdaB, deltaIvalue, itemax, Kextend, bcmax, corenum )
% This function enumerate all the given parameter sets on image series moviethis
%   trjansenum, mvfullenum, trjansdetailenum are cell arrays.
Tmax = size(moviethis,3);
N = size(sprate,1);
trjansenum = cell(N,1);
mvfullenum = cell(N,1);
trjansdetailenum = cell(N,1);

%48 hours for parallel workers.
TimeOut=2880;
parpool('local',corenum,'IdleTimeout',TimeOut);

parfor i=1:N
    lambda = [lambdaA(i) + [0:Tmax-1] * lambdaA(i)*lambdaB(i), 0];
    deltaI = deltaIvalue*ones(1,Tmax);
    [ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), sig, psfdecay, bsize, D, alpha, lambda, deltaI, sprate(i), itemax, Kextend, bcmax, trjtrue, 'no', ['paratunning',num2str(i)] );
    trjansenum{i} = trjans;
    mvfullenum{i} = mvfull;
    trjansdetailenum{i} = trjansdetail;
end

delete(gcp('nocreate'));

end

