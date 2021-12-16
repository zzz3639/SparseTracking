function [ mvgrad ] = grad_t_mv( mvori )
%GRAD_T_MV Summary of this function goes here
%   Detailed explanation goes here

s1 = size(mvori,1);
s2 = size(mvori,2);
T = size(mvori,3);
mvgrad = zeros(s1,s2,T-1);
for i=1:T-1
    mvgrad(:,:,i) = mvori(:,:,i+1)-mvori(:,:,i);
end

end

