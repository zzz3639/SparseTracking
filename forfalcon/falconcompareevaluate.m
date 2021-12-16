function [ Accmatch, Acclink2f, Acclink3f ] = falconcompareevaluate( falconans, trjtruefull )
%FALCONCOMPAREEVALUATE Summary of this function goes here
%   Detailed explanation goes here
N = length(falconans);
oneonematchnum = zeros(N,1);
linkingmatch2f = zeros(N,1);
molnum = zeros(N,1);
linkingnum2f = zeros(N,1);
linkingmatch3f = zeros(N,1);
linkingnum3f = zeros(N,1);
thint = 200;
for i=1:N
    falconansthis = falconans{i};
    trjtruethis = trjtruefull{i};
    trjfalconremain = intensityfilter4column(falconansthis,thint);
    [match1, match2] = oneonematch(trjtruethis.track, trjfalconremain);
    oneonematchnum(i) = sum(match1>0);
    molnum(i) = size(trjtruethis.track,1);
    [ linkingmatchthis ] = linkingmatchtrue( trjtruethis.track, trjfalconremain);
    linkingmatch2f(i) = sum(linkingmatchthis>0);
    linkingnum2f(i) = sum(linkingmatchthis>-1);
    [linkingmatchthis3f] = linkingmatchtrue3frame( trjtruethis.track, trjfalconremain);
    linkingmatch3f(i) = sum(linkingmatchthis3f>0);
    linkingnum3f(i) = sum(linkingmatchthis3f>-1);
end

Accmatch = oneonematchnum./molnum;
Acclink2f = linkingmatch2f./linkingnum2f;
Acclink3f = linkingmatch3f./linkingnum3f;

end

