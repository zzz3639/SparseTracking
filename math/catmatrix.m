function [ Mcat ] = catmatrix( M, I1, I2 )
% cat elements from M indexed by I1,I2, I1,I2 are matrices of the same size.
%   size(Mcat) == size(I1) == size(I2), Mcat(i,j) = M(I1(i,j),I2(i,j));

sizeM = size(M);
m = size(I1,1);
n = size(I1,2);
I1rs = reshape(I1,m*n,1);
I2rs = reshape(I2,m*n,1);
IMreshape = sub2ind(sizeM,I1rs,I2rs);
IM = reshape(IMreshape,m,n);
Mcat = M(IM);
end

