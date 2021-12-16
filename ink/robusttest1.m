function [ mcmv ] = robusttest1( trjtrue )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

startpoint = [20.5,20.5]
S = [42,42];
zoom = 10;
mparticle = 20000;
Tmaxmc = 2000;
pjump = 0.1;
knnjump = 5;
trjbad = trjtrue;
trjbad(:,1:2) = trjbad(:,1:2) + 0.05*randn(size(trjbad,1),2);
zoom = 10;

[ reswalking ] = runwalkingink( trjbad, knnjump, 3, startpoint, pjump, mparticle, Tmaxmc );
[ trjmc ] = reswalking2trj( reswalking );
[ mcmv ] = visualizetrjtrack( trjmc, [42,42], 0, zoom ,10);
m = 4;
for i=1:m
    i
    [ reswalking ] = runwalkingink( trjbad, knnjump, 3, startpoint, pjump, mparticle, Tmaxmc );
    [ trjmc ] = reswalking2trj( reswalking );
    [ mcmvtemp ] = visualizetrjtrack( trjmc, [42,42], 0, zoom ,10);
    mcmv = mcmv + mcmvtemp;
end

end

