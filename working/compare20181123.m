clc;
clear;
thint = 200;
thdist = 0.3;
[ ~, ~, ~, ~, ~, ~, bcmax ] = defaultparaset();
% compute accuracy 1
P1 = '/home/zzz/trackingexperiments/datasetbac/newcomputer/500frame/para_02_02/CompareRho1/resCompareRho1fix/resfixresCompareRho1.mat';
load(P1);
trjans = trjansdetailtotrjans(trjansdetail, bcmax, thint);
[m,s] = linkingmatch( trjtrue.track, trjans.track, thdist );
acc1 = m/s;
Ntrue = size(trjtrue.track,1);
[m1,m2] = oneonematch(trjtrue.track,trjans.track);
accm1 = sum(m1>0)/Ntrue;
[ lm ] = linkingmatchtrue( trjtrue.track, trjans.track);
acclm1 = sum((lm>0))/sum(lm>-1);
% compute accuracy 2
P2 = '/home/zzz/trackingexperiments/datasetbac/newcomputer/500frame/para_02_02_b/CompareRho1/resCompareRho1fix/resfixresCompareRho1.mat';
load(P2);
trjans = trjansdetailtotrjans(trjansdetail, bcmax, thint);
[m,s] = linkingmatch( trjtrue.track, trjans.track, thdist );
acc2 = m/s;
Ntrue = size(trjtrue.track,1);
[m1,m2] = oneonematch(trjtrue.track,trjans.track);
accm2 = sum(m1>0)/Ntrue;
[ lm ] = linkingmatchtrue( trjtrue.track, trjans.track);
acclm2 = sum((lm>0))/sum(lm>-1);
% compute accuracy 3
P3 = '/home/zzz/trackingexperiments/datasetbac/CompareParaFix500frame/CompareRho1fix/resfixCompareRho1.mat';
load(P3);
trjans = trjansdetailtotrjans(trjansdetail, bcmax, thint);
[m,s] = linkingmatch( trjtrue.track, trjans.track, thdist );
acc3 = m/s;
Ntrue = size(trjtrue.track,1);
[m1,m2] = oneonematch(trjtrue.track,trjans.track);
accm3 = sum(m1>0)/Ntrue;
[ lm ] = linkingmatchtrue( trjtrue.track, trjans.track);
acclm3 = sum((lm>0))/sum(lm>-1);