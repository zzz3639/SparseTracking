clear;
clc;
P = '/home/zzz/trackingexperiments/datasetbac/';
Nset = 8;
Im = 5;
thint = 200;
thdist = 0.3;
trjMTTcell = cell(Nset,1); %5 columns
trjtruecell = cell(Nset,1);
AccMTT = zeros(Nset,1);
AcclinkingMTT = zeros(Nset,1);
RecallNum = zeros(Nset,1);
for i=1:8
    load([P,'CompareRho',num2str(i),'.mat']);
    [trjansMTT, linkingsMTT] = runMTTtracking(moviethis);
    [ trjMTT ] = MTTlinks( linkingsMTT );
    trjMTT(:,3) = trjMTT(:,3)*Im;
    trjMTTfilted = intensityfilter5column(trjMTT,thint);
    trjMTTcell{i} = trjMTTfilted;
    trjtruecell{i} = trjtrue;
    [m1,m2] = oneonematch( trjtrue.track, trjMTTfilted );
    a = length(m1); b = sum(m1>0);
    AccMTT(i) = b/a;
    [a,b] = linkingmatch( trjtrue.track, trjMTTfilted, thdist );
    AcclinkingMTT(i) = a/b;
    RecallNum(i) = b;
end 

save MTTresults.mat;
