function [ trjfilted ] = intensityfilter4column( trj, threshold )
%intensityfilter4column removes molecules whose intensity below threshold
%  inputs:
%   trj: [x,y,I,t], before filter
%   threshold: non-negative number
%  outputs:
%   trjfilted: [x,y,I,t], after filter

v = find(trj(:,3)>threshold);
trjfilted = trj(v,:);

end

