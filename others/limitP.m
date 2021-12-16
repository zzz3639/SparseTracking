function [ r, Mwrong ] = limitP( mu, sigma, boardsize, Tmax )
%LIMITP Summary of this function goes here
% this is a ring world
%   mu, molecule density
%   sigma: molecule moving velocity, gaussian sigma
%   boardsize: the board doing simulation
%   Tmax: How long the simulatio is going.

rng('shuffle');
N = ceil(boardsize*boardsize*mu);

% initialize N particles at random %
M = rand(N,2)*boardsize;

% run the brownian motion engine
Mt = zeros(N,2,Tmax+1);
Mt(:,:,1) = M;
for i=1:Tmax
    Mt(:,:,i+1) = Mt(:,:,i) + randn(N,2)*sigma;
    Mt(:,:,i+1) = Mt(:,:,i+1) - (Mt(:,:,i+1)>boardsize)*boardsize;
end

Mwrong = zeros(N,Tmax);

for i=1:Tmax
    D = zeros(N,N);
    Dx = zeros(N,N,3);
    Dy = zeros(N,N,3);
    Dx(:,:,1) = (repmat(Mt(:,1,i),1,N) - repmat([Mt(:,1,i+1)]',N,1)).^2;
    Dx(:,:,2) = (repmat(Mt(:,1,i),1,N) - repmat([Mt(:,1,i+1)]',N,1) - boardsize).^2;
    Dx(:,:,3) = (repmat(Mt(:,1,i),1,N) - repmat([Mt(:,1,i+1)]',N,1) + boardsize).^2;
    Dy(:,:,1) = (repmat(Mt(:,2,i),1,N) - repmat([Mt(:,2,i+1)]',N,1)).^2;
    Dy(:,:,2) = (repmat(Mt(:,2,i),1,N) - repmat([Mt(:,2,i+1)]',N,1) - boardsize).^2;
    Dy(:,:,3) = (repmat(Mt(:,2,i),1,N) - repmat([Mt(:,2,i+1)]',N,1) + boardsize).^2;
    Dxnew = min(Dx,[],3);
    Dynew = min(Dy,[],3);
    D = Dxnew+Dynew;
    Dfalse = D + D';
    Dtrue = repmat(diag(D),1,N);
    Dtrue = Dtrue + Dtrue';
    Dthis = Dtrue - Dfalse;
    [u,v] = max(Dthis,[],2);
    Mwrong(find(u>0),i) = 1;
end

r = sum(sum(Mwrong));
r = r/(N*Tmax);

end

