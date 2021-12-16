function [ idxsub ] = subsetsearch( vfull, vsub )
% vfull and vsub are vector of positive integers
% vsub is a subset of vfull, return the index of elements of vsub in vfull

nsub = length(vsub);
idxsub = zeros(nsub,1);
for i=1:nsub
    idxsub(i) = find(vfull==vsub(i));
end

end

