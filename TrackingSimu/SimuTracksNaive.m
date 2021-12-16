function [Movies, Trj] = SimuTracksNaive(Smol,S,I0,DI0,alpha,Dmu,Rho,Tmax,Nrepeat,Mrepeat,no,sigma,obnoise,pos0)
%Modified in 2017.11.06 by ZHANG Haowen
%Usage: [Movies, Trj] = SimuTracksNaive(Smol,S,I0,DI0,alpha,Dmu,Rho,Tmax,Nrepeat,Mrepeat,no,sigma,obnoise)
%  Smol: [smol1,smol2,d1,d2], molecule activated inside area (d1:d1+smol1-1,d2:d2+smol2-1)
%  S: [s1,s2], area of focus at imaging
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alpha: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: diffusion parameter, normal, sigma/pixels
%  Rho: The distribution of activated molecules in frame 0, a matrix of size [smol1,smol2];
%  Tmax: total number of frames simulated, first frame 0, last frame Tmax-1
%  Nrepeat: simulate Nrepeat tracks
%  Mrepeat: Mrepeat imagings for one trajectory
%  no: noise, count by photon number
%  sigma: PSF width
%  obnoise: Additional noise, for each pixel, signal=photoncount+sqrt(photoncount)*randn()*obnoise
%  pos0: molecule positions in the first frame, if ~exist('pos0'), generate at random


rng('shuffle');
% define image positions list
s1 = Smol(1);
s2 = Smol(2);
cx1 = [0:s1-1]';
cx1 = repmat(cx1,1,s2);
cx2 = [0:s2-1];
cx2 = repmat(cx2,s1,1);
cx = [reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];
LogI0 = log(I0);
%initiate Pic0
if exist('pos0')
    NAC = size(pos0,1);
    Pic0 = zeros(NAC,3);
    Pic0(:,1:2) = pos0;
    Pic0(:,3) = exp(LogI0+DI0*randn(NAC,1));    
else
    MAC = poissrnd(Rho);
    NAC = sum(sum(MAC));
    Pic0 = zeros(NAC,3);
    Pic0(:,3) = exp(LogI0+DI0*randn(NAC,1));
    MACline = reshape(MAC,s1*s2,1);
    V = find(MACline);
    j=0;
    for i=1:length(V)
        MACline(V(i));
        Pic0(j+1:j+MACline(V(i)),1:2) = repmat(cx(V(i),:),MACline(V(i)),1) + rand(MACline(V(i)),2) - 0.5;
        j = j+MACline(V(i));
    end
    Pic0(:,1:2) = Pic0(:,1:2) + repmat([Smol(3),Smol(4)],NAC,1);
end

Movies = cell(Nrepeat,Mrepeat);

Trj = SiliconMovementNaive(I0,DI0,alpha,Dmu,Smol(1:2),Pic0,Tmax,Nrepeat);
for i=1:Nrepeat
    trjthis = Trj{i};
    Moviethis = SiliconImagingTracks( Mrepeat, S, trjthis, no, sigma, obnoise );
    Movies(i,:) = Moviethis';
end

end
