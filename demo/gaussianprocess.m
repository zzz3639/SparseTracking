function [ X ] = gaussianprocess( N, T, alpha )
%GAUSSIANPROCESS Summary of this function goes here
%   X: (T+1)*N
X = randn(1,N);

for i=1:T
    Xnext = alpha*X(end,:) + sqrt(1-alpha*alpha)*randn(1,N);
    X = [X;Xnext];
end

end

