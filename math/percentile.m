function [ P ] = percentile( A, percent )
% 75% percentile of an array
%   Detailed explanation goes here

    A = sort(A);
    N = ceil(percent*length(A));
    P = A(N);
end

