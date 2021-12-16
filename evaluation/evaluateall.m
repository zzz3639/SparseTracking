function [ Acc, Acclinktrue2f, Acclinktrue3f, Acclink2f, Acclink3f, linkrecall ] = evaluateall( trjtrue, trjans, intensityth, thdist )
%EVALUATEALL Summary of this function goes here
% outputs:
%   Acc: Detection accuracy
%   Acclinktrue2f: linking accuracy estimated solely by localiztions
%   Acclinktrue3f: linking accuracy estimated solely by localiztions, 3 frames
%   Acclink: linking accuracy
%   Acclink3f: linking accuracy, 3 frames

matchdata = zeros(1,3);
linkingmatchnum = zeros(1,2);
linkingmatchnum3f = zeros(1,2);
linking = zeros(1,2);
if size(trjans,2)==5
    trjansres = intensityfilter5column(trjans,intensityth);
else
    trjansres = intensityfilter4column(trjans,intensityth);
end
[match1,match2] = oneonematch( trjtrue, trjansres );
matchdata(1) = sum(match1>0);
matchdata(2) = length(match1);
matchdata(3) = length(match2);
[ linkingmatch2f ] = linkingmatchtrue( trjtrue, trjansres);
linkingmatchnum(1) = sum(linkingmatch2f>0);
linkingmatchnum(2) = sum(linkingmatch2f>-1);
[ linkingmatch3f ] = linkingmatchtrue3frame( trjtrue, trjansres);
linkingmatchnum3f(1) = sum(linkingmatch3f>0);
linkingmatchnum3f(2) = sum(linkingmatch3f>-1);

Acc = matchdata(1)./matchdata(2);
Acclinktrue2f = linkingmatchnum(1)./linkingmatchnum(2);
Acclinktrue3f = linkingmatchnum3f(1)./linkingmatchnum3f(2);

linkrecall = zeros(2,1);
if size(trjansres,2)==5
    [ m,s ] = linkingmatch( trjtrue, trjansres, thdist);
    linking(1) = m;
    linking(2) = s;
    Acclink2f = linking(1)./linking(2);
    linkrecall(1) = s;
    [ m,s ] = linkingmatchkframe( trjtrue, trjansres, thdist, 3);
    linking(1) = m;
    linking(2) = s;
    Acclink3f = linking(1)./linking(2);
    linkrecall(2) = s;
else
    Acclink2f = 0;
    Acclink3f = 0;
    linkrecall(:,1) = 0;
end

end

