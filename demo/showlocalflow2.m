function [ vximg, vyimg, numimg ] = showlocalflow2( trj, focusarea, zoom, scalefactor )
% given a trajectory matrix [x,y,I,t,id], show local diffusion rate of each grid in focus area
% usage:
%   inputs:
%     trj: [x,y,I,t,id]
%     focusarea: [xs,xe,ys,ye], size of a grid in this area is 1*1
%     zoom: size of grid shown in output diffimg.
%     scalefactor: to plot opticalflow, larger value gives longer arrow,
%         give negative number if you don't want to plot
%   outputs:
%     vximg,vyimg: size (xs-xe)*(ys-ye), x,y velocity
%     numimg: number of links in each grid

xs = focusarea(1);  xe = focusarea(2);
ys = focusarea(3);  ye = focusarea(4);
trjzoom = trj;
L = size(trj,1);
trjzoom(:,1:2) = (trj(:,1:2)-repmat([xs,ys],L,1))*zoom;

xezoom = floor((xe-xs)*zoom);
yezoom = floor((ye-ys)*zoom);
[ vximg, vyimg, numimg ] = showlocalflow( trjzoom, [0,xezoom,0,yezoom] );
vximg = vximg/zoom;
vyimg = vyimg/zoom;

if scalefactor>0
    flow = opticalFlow(vyimg,vximg);
    plot(flow, 'ScaleFactor', scalefactor);
end

end

