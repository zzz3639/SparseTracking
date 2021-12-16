function [Movies, Trj] = SimuTracksPoisson(Smol,S,AC,SV,I0,DI0,alpha,Dmu,Rho,Pic0,Tmax,Nrepeat,Mrepeat,no,sigma,obnoise)
%Modified in 2018.12.04 by ZHANG Haowen
% add poisson noise(Inew=poissrnd(I)) to emitter intensities.
%Usage: [Movies] = SimuTracksPoisson(Smol,S,AC,SV,I0,DI0,alpha,Dmu,Rho,Pic0,Tmax,Nrepeat,Mrepeat,no,sigma,obnoise)
% position of first pixel is (0,0)
% inputs:
%  Smol: [smol1,smol2,d1,d2], molecule activated inside area (d1:d1+smol1-1,d2:d2+smol2-1)
%  S: [s1,s2], area of focus at imaging
%  AC: activation rate, AC*Rho_ij is the number of molecules activated of pixel area ij per frame
%  SV: average survival length of each activated molecule
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alpha: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: diffusion parameter, normal, sigma/pixels
%  Rho: The distribution of activated molecules, is normalized inside this function
%      matrix of size s1*s2
%  Pic0: activated molecules in frame 0, [X,Y,I]
%  Tmax: total number of frames simulated, first frame 0, last frame Tmax-1
%  Nrepeat: simulate Nrepeat tracks
%  Mrepeat: Mrepeat imagings for one trajectory
%  no: noise, count by photon number
%  sigma: PSF width
%  obnoise: Additional noise, for each pixel, signal=photoncount+sqrt(photoncount)*randn()*obnoise
% outputs:
%  Movies:
%  Trj:

Movies = cell(Nrepeat,Mrepeat);
Rho = Rho/sum(sum(Rho));
Trj = SiliconMovement(AC,SV,I0,DI0,alpha,Dmu,Smol(1:2),Rho,Pic0,Tmax,Nrepeat);
for i=1:Nrepeat
    trjthis = Trj{i};
    trjthis(:,3) = poissrnd(trjthis(:,3));
    Trj{i} = trjthis;
end
for i=1:Nrepeat
    trjthis = Trj{i};
    trjthis(:,1:2) = trjthis(:,1:2) + [Smol(3)*ones(size(trjthis,1),1),Smol(4)*ones(size(trjthis,1),1)];
    Trj{i} = trjthis;
    Moviethis = SiliconImagingTracks( Mrepeat, S, trjthis, no, sigma, obnoise );
    Movies(i,:) = Moviethis';
end

end
