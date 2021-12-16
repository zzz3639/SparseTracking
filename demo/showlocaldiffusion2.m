function [ diffimg, numimg ] = showlocaldiffusion2( trj, zoom, focusarea )
% given a trajectory matrix [x,y,I,t,id], show local diffusion rate of each grid in focus area
% usage:
%   inputs:
%     trj: [x,y,I,t,id]
%     focusarea: [xs,xe,ys,ye], area to be plotted, same coordinate with
%         x,y in trj
%     zoom: [x,y]=[x,y]*zoom, during statistic
%   outputs:
%     diffimg: size (xs-xe)zoom*(ys-ye)zoom, value equals to msd of each grid
%     numimg: number of links in each grid

xs = focusarea(1);
xe = focusarea(2);
ys = focusarea(3);
ye = focusarea(4);
trj(:,1) = (trj(:,1)-xs)*zoom;
trj(:,2) = (trj(:,2)-ys)*zoom;
R = [[xe-xs],[ye-ys]]*zoom;
R = floor(R);
[ diffimg, numimg ] = showlocaldiffusion( trj, [0,R(1),0,R(2)] );
diffimg = diffimg/zoom/zoom;

end