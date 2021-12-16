function [ diffimg, numimg ] = showlocaldiffusion( trj, focusarea )
% given a trajectory matrix [x,y,I,t,id], show local diffusion rate of each grid in focus area
% usage:
%   inputs:
%     trj: [x,y,I,t,id]
%     focusarea: [xs,xe,ys,ye], size of a grid in this area is 1*1
%     zoom: size of grid shown in output diffimg.
%   outputs:
%     diffimg: size (xs-xe)zoom*(ys-ye)zoom, value equals to msd of each grid
%     numimg: number of links in each grid

[trjunique,lbnum]= molidunique(trj);
trjsort = sorttrackid(trjunique);
xs = focusarea(1);
xe = focusarea(2);
ys = focusarea(3);
ye = focusarea(4);
lx = floor(xe-xs);
ly = floor(ye-ys);
diffimg = zeros(lx,ly);
numimg = zeros(lx,ly);
L = zeros(lbnum,1);
for i=0:lbnum-1
    L(i+1) = length(find(trjsort(:,5)==i));
end

k=1;
for i=1:lbnum
    for j=1:L(i)-1
        gx = ceil(trjsort(k,1)-xs);
        gy = ceil(trjsort(k,2)-ys);
        if gx>0&&gx<lx+1&&gy>0&&gy<ly+1
            numimg(gx,gy) = numimg(gx,gy) + 1;
            sdthis = (trjsort(k+1,1)-trjsort(k,1))^2 + (trjsort(k+1,2)-trjsort(k,2))^2;
            diffimg(gx,gy) = diffimg(gx,gy) + sdthis;
        end
        k=k+1;
    end
    k=k+1;
end

for i=1:lx
    for j=1:ly
        if numimg(i,j)==0
        else
            diffimg(i,j) = diffimg(i,j)/numimg(i,j);
        end
    end
end

end

