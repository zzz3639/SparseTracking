function [ fknzarray ] = firstknonzero( M, k )
% find first k non-zero element for every line of M
%   M: matrix, m*n
%   fknzarray: m*k, index
S = size(M);
fknzarray = zeros(S(1),k);
for i=1:k
    fnzarray = firstnonzero(M);
    fknzarray(:,i) = fnzarray;
    %add one column to M
    Mb = [M,zeros(S(1),1)];
    Mb(sub2ind(size(Mb), [1:S(1)], fnzarray')) = 0;
    %change M
    M = Mb(:,1:S(2));
end

end

