clear;
clc;
Nset=8;
intensityth=200;
P = '/home/zzz/trackingexperiments/datasetbac/newserver/500frame/para_02_02/'
matchdata = zeros(Nset,3);
linkingmatchnum = zeros(Nset,2);
linkingmatchnum3f = zeros(Nset,2);
%load match numbers and save the numbers to matches
for i=1:Nset
    load([P,'CompareRho',num2str(i),'/resCompareRho',num2str(i),'fix/resfixresCompareRho',num2str(i),'.mat']);
    trjansres = intensityfilter4column(trjans.track,intensityth);
    [match1,match2] = oneonematch( trjtrue.track, trjansres );
    matchdata(i,1) = sum(match1>0);
    matchdata(i,2) = length(match1);
    matchdata(i,3) = length(match2);
    [ linkingmatch2f ] = linkingmatchtrue( trjtrue.track, trjansres);
    linkingmatchnum(i,1) = sum(linkingmatch2f>0);
    linkingmatchnum(i,2) = sum(linkingmatch2f>-1);
    [ linkingmatch3f ] = linkingmatchtrue3frame( trjtrue.track, trjansres);
    linkingmatchnum3f(i,1) = sum(linkingmatch3f>0);
    linkingmatchnum3f(i,2) = sum(linkingmatch3f>-1);
end
Acc = zeros(Nset,1);
Acclink2f = zeros(Nset,1);
Acclink3f = zeros(Nset,1);
Acc = matchdata(:,1)./matchdata(:,2);
Acclink2f = linkingmatchnum(:,1)./linkingmatchnum(:,2);
Acclink3f = linkingmatchnum3f(:,1)./linkingmatchnum3f(:,2);
