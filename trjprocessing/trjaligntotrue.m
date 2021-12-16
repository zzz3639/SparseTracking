function [ linealigned ] = trjaligntotrue( trjtrue, trjans )
% align every emitter in trjans to trjtrue;
% inputs: trjtrue: [x,y,I,t,id]
%         trjans: at least 4 columns: [x,y,I,t]
% outputs: 
%    linealigned: same line number with trjans, i to linealigned(i) in trjtrue

T = max(trjans(:,4));
N = size(trjans,1);
linealigned = zeros(N,1);

for t=0:T
    vans = find(trjans(:,4)==t);
    vtrue = find(trjtrue(:,4)==t);
    trjanstime = trjans(vans,:);
    trjtruetime = trjtrue(vtrue,:);
    [ framealigned ] = MolListAlign(trjanstime(:,1:3), trjtruetime(:,1:3));
    linealigned(vans,1) = vtrue(framealigned(:,7),1);
end

end

