function [ spratebest, lambdaAbest, lambdaBbest, matchnum, trjansenum, trjans5cenum ] = runparatunning( moviethis, trjtrue, sig, psfdecay, bsize, D, alpha, thint, corenum )
%RUNPARATUNNING Summary of this function goes here
%   thint: dim particles under thint are removed
%   corenum: number of computer cores for parfor

[ sprate, lambdaA, lambdaB, deltaIvalue, itemax, Kextend, bcmax ] = defaultparaset();
[ trjansenum, ~, trjansdetailenum ] = enumparaset( moviethis, trjtrue, sig, psfdecay, bsize, D, alpha, sprate, lambdaA, lambdaB, deltaIvalue, itemax, Kextend, bcmax, corenum );

[ matchnum ] = enumacc( trjtrue, trjansenum, thint );

Mpara = length(sprate);

trjans5cenum = cell(Mpara,1);

for i=1:Mpara
    trjans5cenum{i} = trjansdetailtotrjans(trjansdetailenum{i},bcmax,thint);
end

matches = zeros(Mpara,3);

matches(:,1) = matchnum;
matches(:,2) = size(trjtrue.track,1);
for i=1:Mpara
    matches(i,3) = size(trjansenum{i}.track,1);
end

Mnum = matches(:,1)./(matches(:,2)+matches(:,3));

[u,v] = max(Mnum);

spratebest = sprate(v);
lambdaAbest = lambdaA(v);
lambdaBbest = lambdaB(v);

end

