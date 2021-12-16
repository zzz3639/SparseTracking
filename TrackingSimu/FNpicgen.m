function [ pos ] = FNpicgen( Smol, N, Rho )
%Fixed number particle positions generater
%Usage:  pos = FNpicge( Smol, N, Rho )
%  Smol: [smol1,smol2,d1,d2], molecule activated inside area (d1:d1+smol1-1,d2:d2+smol2-1)
%  N: number of molecules to be placed
%  Rho: The distribution of activated molecules in frame 0, a matrix of size [smol1,smol2];

% define image positions list
rng('shuffle');
s1 = Smol(1);
s2 = Smol(2);
cx1 = [0:s1-1]';
cx1 = repmat(cx1,1,s2);
cx2 = [0:s2-1];
cx2 = repmat(cx2,s1,1);
cx = [reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];

%compute pos
Rho = Rho/sum(sum(Rho));
k=0;
pos = zeros(0,2);
while k<N
    MAC = poissrnd(N*Rho);
    NAC = sum(sum(MAC));
    Pic0 = zeros(NAC,2);
    MACline = reshape(MAC,s1*s2,1);
    V = find(MACline);
    j=0;
    for i=1:length(V)
        Pic0(j+1:j+MACline(V(i)),1:2) = repmat(cx(V(i),:),MACline(V(i)),1) + rand(MACline(V(i)),2) - 0.5;
        j = j+MACline(V(i));
    end
    pos = [pos;Pic0];
    k = k+NAC;
end
M = size(pos,1);
Mperm = randperm(M);
pos = pos(Mperm(1:N),:);
pos = pos + repmat([Smol(3),Smol(4)],N,1);

end

