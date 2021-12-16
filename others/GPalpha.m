function [ alpha ] = GPalpha( DI0, rho )
%GPALPHA Summary of this function goes here
%   Detailed explanation goes here
rng('shuffle');
T = 10000;
M = zeros(T,1);
M(1) = DI0*randn();
for i=2:T
    M(i) = rho*M(i-1)+sqrt(1-rho*rho)*randn()*DI0;
end
alpha = mean(abs(M(2:T,1)-M(1:T-1,1)));
end

