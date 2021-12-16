function [ trj ] = reswalking2trj( reswalking )
%RESWALKING2TRJ Summary of this function goes here
%   Detailed explanation goes here

m = size(reswalking,1);
Tmax = size(reswalking,3);

trj = zeros(m*Tmax,5);
for i=1:Tmax
    trj((i-1)*m+1:i*m,1:2) = reswalking(:,:,i);
    trj((i-1)*m+1:i*m,3) = ones(m,1);
    trj((i-1)*m+1:i*m,4) = i-1;
    trj((i-1)*m+1:i*m,5) = [0:m-1]';
end

end

