function [mvplay] = trjplay3( trjtrue, trjans, trjans2, S, zoom, brightness )
%A movie player shows the the trjtrue
%  inputs:
%   trjtrue: at least 4 columns [x,y,I,t]
%   trjans: at least 4 columns [x,y,I,t]
%   S: [s1,s2]
%   zoom: zoom rate, movie size equals to S*zoom
%   brightness: 1 by default

Tmax1 = max(trjtrue(:,4));
Tmax2 = max(trjans(:,4));
Tmax3 = max(trjans2(:,4));
Tmax = max([Tmax1,Tmax2,Tmax3])+1;
mvtrue = visualizetrjtrack(trjtrue,S,Tmax,zoom,1);
mvans = visualizetrjtrack(trjans,S,Tmax,zoom,1);
mvans2 = visualizetrjtrack(trjans2,S,Tmax,zoom,1);

m1 = max(max(max(mvtrue)));
m2 = max(max(max(mvans)));
m3 = max(max(max(mvans2)));
m = max([m1,m2,m3]);

s1 = S(1);
s2 = S(2);
mvplay = zeros(s1*zoom,3*s2*zoom,Tmax);
for i=1:Tmax
    mvplay(:,:,i) = [mvtrue(:,:,i),mvans(:,:,i),mvans2(:,:,i)];
end
clear mvtrue mvans mvans2
mvplay = mvplay*brightness/m;
implay(mvplay);


end

