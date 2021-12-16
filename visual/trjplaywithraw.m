function [mvplay] = trjplaywithraw( moviethis, trjtrue, trjans, S, zoom, brightness )
%A movie player shows the the trjtrue
%  inputs:
%   moviethis: raw movie
%   trjtrue: at least 4 columns [x,y,I,t]
%   trjans: at least 4 columns [x,y,I,t]
%   S: [s1,s2]
%   zoom: zoom rate, movie size equals to S*zoom
%   brightness: 1 by default

Tmaxmv = size(moviethis,3);
Tmax1 = max(trjtrue(:,4));
if length(Tmax1)==0
    Tmax1 = 0;
end
Tmax2 = max(trjans(:,4));
if length(Tmax2)==0
    Tmax2 = 0;
end
Tmax = max(Tmax1,Tmax2)+1;
Tmax = min(Tmax,Tmaxmv);
mvtrue = visualizetrjtrack(trjtrue,S,Tmax,zoom,1);
mvans = visualizetrjtrack(trjans,S,Tmax,zoom,1);

m1 = max(max(max(mvtrue)));
m2 = max(max(max(mvans)));
m = max(m1,m2);

mvplayans = cat(2,mvtrue,mvans)*brightness/m;
clear mvtrue;
clear mvans;
mvplayraw = zeros([S*zoom,Tmax]);
mraw = max(max(max(moviethis)));
for i=1:Tmax
    mvplayraw(:,:,i) = kron([moviethis(:,:,i)]',ones(zoom,zoom))/mraw;
end
mvplay = cat(2,mvplayraw,mvplayans);
clear mvplayraw;
clear mvplayans;
implay(mvplay);


end

