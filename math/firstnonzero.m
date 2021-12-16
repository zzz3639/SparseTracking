function [ fnzarray ] = firstnonzero( M )
% find first non-zero element for every line of M
%   M: matrix, m*n
%   fnzarray: m*1, index
S = cumsum(M,2);
fnzarray = sum((S==0),2)+1;

end

