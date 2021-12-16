clear;
clc;
Nset=8;
Mpara=6;
intensityth=200;
P = '/home/zzz/trackingexperiments/datasetbac/Compare03ParaSelect100frame/'
matches = zeros(Mpara,Nset,3);
linkingmatches = zeros(Mpara,Nset,2);
bestpara = zeros(Nset,3);
%load match numbers and save the numbers to matches
for i=1:Nset
    load([P,'Compare03Rho',num2str(i),'/resCompare03Rho',num2str(i),'.mat']);
    matches(:,i,1) = matchnum;
    matches(:,i,2) = size(trjtruetunning.track,1);
    for j=1:Mpara
        matches(j,i,3) = size(trjansenum{j}.track,1);
        trackfilted = intensityfilter4column( trjansenum{j}.track, intensityth );
        linkingmatchthis = linkingmatchtrue( trjtruetunning.track, trackfilted);
        linkingmatches(j,i,1) = sum(linkingmatchthis>0);
        linkingmatches(j,i,2) = sum(linkingmatchthis>-1);
    end
end
Mnum = matches(:,:,1)./(matches(:,:,2)+matches(:,:,3));
[ sprate, lambdaA, lambdaB ] = defaultparaset();
Acc = zeros(Nset,1);
Acclink = zeros(Nset,3);
for i=1:Nset
    Mnumarray = Mnum(:,i);
    [u,v] = max(Mnumarray);
    bestpara(i,1) = sprate(v);
    bestpara(i,2) = lambdaA(v);
    bestpara(i,3) = lambdaB(v);
    Acc(i) = matches(v,i,1)/matches(v,i,2);
    Acclink(i,1) = linkingmatches(v,i,1)/linkingmatches(v,i,2);
    Acclink(i,2) = linkingmatches(v,i,1);
    Acclink(i,3) = linkingmatches(v,i,2);
end

