function [ alpha ] = alphaest( trackin )
%ALPHAEST Summary of this function goes here
%   Detailed explanation goes here

trackin = sorttrackid(trackin)

alpha=0;
N=0;
L=size(trackin,1);

for i=2:L
    if trackin(i,5)==trackin(i-1,5)
        N = N+1;
        alpha = alpha + abs(log(trackin(i,3))-log(trackin(i-1,3)));
    end
end

alpha = alpha/N;

end

