function [ trackR, trackP ] = trackingerror1( trjans , trjtrue, IntFilter )
%count recall and precision between trjans and trjtrue, 
%only molecules whose intensity larger than IntFilter are counted
%  Usage: [ trackR, trackP ] = trackingerror1( trjans , trjtrue, IntFilter )
trjanstrack = trjans.track;
trjtruetrack = trjtrue.track;
v1 = find(trjanstrack(:,3)>IntFilter);
v2 = find(trjtruetrack(:,3)>IntFilter);
trjanstrack = trjanstrack(v1,:);
trjanstrack(:,1:2) = trjanstrack(:,[2,1]);
trjtruetrack = trjtruetrack(v2,:);
Tmax = max(max(trjanstrack(:,4)),max(trjanstrack(:,4)))+1;
trackR = zeros(Tmax,3);
trackP = zeros(Tmax,3);
for t=1:Tmax
    v1 = find(trjanstrack(:,4)==t-1);
    v2 = find(trjtruetrack(:,4)==t-1);
    trjansthis = trjanstrack(v1,:);
    trjtruethis = trjtruetrack(v2,:);
    algn1 = MolListAlign(trjansthis, trjtruethis);
    algn2 = MolListAlign(trjtruethis, trjansthis);
    trackR(t,1) = size(trjtruethis,1);
    trackP(t,1) = size(trjansthis,1);
    trackR(t,2) = length(unique(algn1(:,end)));
    trackP(t,2) = length(unique(algn2(:,end)));
    trackR(t,3) = mean(algn1(:,4));
    trackP(t,3) = mean(algn2(:,4));
end

end

