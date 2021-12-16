function [ trjMTT ] = MTTlinks( linkings )
%change MTT trajectory into trajectory matrix [x,y,I,t,id]
%   linkings, a cell array of MTT tab_param, merge them into a single trajectory matrix

Tmax = length(linkings) + 1;
frame1 = linkings{1}(:,129:131);
v = find(frame1(:,3)>0);
frame1 = frame1(v,:);
m = length(v);
frame1 = [frame1,zeros(size(frame1,1),1),[0:m-1]'];
trjMTT = frame1;

for t=2:Tmax
    lktemp = linkings{t-1};
    el1 = lktemp(:,129:131);
    el2 = lktemp(:,136:138);
    v = find(el2(:,3)>0);
    eltp1add = el2(v,:);
    eltadd = el1(v,:);
    vl = find(eltadd(:,3)>0);
    eltfrom = eltadd(vl,:);
    vt = find(trjMTT(:,4)==t-2);
    eltto = trjMTT(vt,:);
    [ alignfrom, alignto ] = emitteralias(eltfrom(:,1:3), eltto(:,1:3), 1, 1000, 1);
    link = zeros(size(eltp1add,1),1);
    link(vl(alignfrom)) = vt(alignto);
    framenext = [eltp1add,(t-1)*ones(size(eltp1add,1),1),[1:size(eltp1add,1)]',link];
    trjMTT = linktwotrj( trjMTT, framenext );
end

trjMTT = trjMTT(:,[2,1,3,4,5]);

end

