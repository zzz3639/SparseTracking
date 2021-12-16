function [ mvfull ] = showtrjans( S, trjtrue, trjans )
%Usage: [ mvfull ] = showtrjans( S, trjtrue, trjans )
%   S = [s1,s2]
%   trjtrue: true answer
%   trjans: solver returned answer
zoom = 20;

s1 = S(1);

Ithis = mean(trjtrue.track(:,3));

Tmax = max(trjtrue.track(:,4))+1;

[ vtrue ] = visualizetrjtrack( [trjtrue.track(:,2:-1:1),trjtrue.track(:,3)/Ithis,trjtrue.track(:,4)], s1, zoom, 1);

[ vans ] = visualizetrjtrack( [trjans.track(:,1:2),trjans.track(:,3)/Ithis,trjans.track(:,4)], s1, zoom, 1);

lmv = s1*zoom;

mvfull = zeros(lmv,lmv*2,Tmax);
for i=1:Tmax
    mvfull(1:lmv,1:lmv,i) = vtrue(:,:,i);
    mvfull(1:lmv,lmv+1:lmv*2,i) = vans(:,:,i);
end

end

