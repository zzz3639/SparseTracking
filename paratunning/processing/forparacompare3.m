clear;
clc;
Nset=8;
intensityth=200;
P = '/home/zzz/trackingexperiments/datasetbac/CompareParaFix500frame/'
matchnum = zeros(Nset,3);
linkingmatchnum = zeros(Nset,2);
Acc = zeros(Nset,1);
Acclink = zeros(Nset,1);
%load match numbers and save the numbers to matches
for i=1:Nset
    load([P,'CompareRho',num2str(i),'fix','/resfixCompareRho',num2str(i),'.mat']);
    trjansres = intensityfilter4column(trjans.track,intensityth);
    [match1,match2] = oneonematch( trjtrue.track, trjansres );
    matchnum(i,1) = sum(match1>0);
    matchnum(i,2) = length(match1);
    matchnum(i,3) = length(match2);
    [ linkingmatch ] = linkingmatchtrue( trjtrue.track, trjansres);
    linkingmatchnum(i,1) = sum(linkingmatch>0);
    linkingmatchnum(i,2) = sum(linkingmatch>-1);
end
Acc = zeros(Nset,1);
Acclink = zeros(Nset,1);
Acc = matchnum(:,1)./matchnum(:,2);
Acclink = linkingmatchnum(:,1)./linkingmatchnum(:,2);
