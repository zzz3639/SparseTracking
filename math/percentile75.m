function [ P75 ] = percentile75( A )
% 75% percentile of an array
%   Detailed explanation goes here

    A = sort(A);
    N = ceil(0.75*length(A));
    P75 = A(N);
end

