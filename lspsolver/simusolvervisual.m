function [moviethis, trjtrue, trjans, mvfull] = simusolvervisual(Smol,S,I0,DI0,alphaGP,Dmu,Rho,Tmax,no,sigma,obnoise, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, mv0, pos0)
%Usage: [moviethis, trjtrue, trjans, mvfull] = simusolvervisual(Smol,S,I0,DI0,alphaGP,Dmu,Rho,Tmax,no,sigma,obnoise, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, mv0, pos0)
%  Smol: [smol1,smol2,d1,d2], molecule activated inside area (d1:d1+smol1-1,d2:d2+smol2-1)
%  S: [s1,s2], area of focus at imaging
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alphaGP: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: diffusion parameter, normal, sigma/pixels
%  Rho: The distribution of activated molecules in frame 0, a matrix of size [smol1,smol2];
%  Tmax: total number of frames simulated, first frame 0, last frame Tmax-1
%  Nrepeat: simulate Nrepeat tracks
%  Mrepeat: Mrepeat imagings for one trajectory
%  no: noise, count by photon number
%  sigma: PSF width
%  obnoise: Additional noise, for each pixel, signal=photoncount+sqrt(photoncount)*randn()*obnoise
%  psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%  bsize: boundary size
%  D: diffusion rate
%  alpha: intensity linkage
%  lambda: LSP norm parameter
%  deltaI: LSP norm parameter
%  itemax: maximal iteration number
%  mv0: mv0.track and mv0.no, if not given or assign to be [], randomly initiated inside this function
%  pos0: molecule positions in the first frame

zoom = 20;
if exist('pos0')
    [Movies, Trj] = SimuTracksNaive(Smol,S,I0,DI0,alphaGP,Dmu,Rho,Tmax,1,1,no,sigma,obnoise,pos0);
else
    [Movies, Trj] = SimuTracksNaive(Smol,S,I0,DI0,alphaGP,Dmu,Rho,Tmax,1,1,no,sigma,obnoise);
end
moviethis = Movies{1,1};
trjtrue.track = Trj{1,1};
trjtrue.no = no*ones(Tmax,1);
smol1 = Smol(1); smol2 = Smol(2);
sarea = smol1*smol2;
if exist('mv0')&&~isempty(mv0)
else
    mv0no = sum(sum(moviethis));
    mv0.no = zeros(Tmax,1);
    mv0.no(:,1) = mv0no(1,1,:);
    Rhomv0 = ones(smol1,smol2)/(smol1*smol2);
    Rhomv0 = ceil(sarea/2)*Rhomv0;
    [Mtemp, Trjtemp] = SimuTracksNaive(Smol,S,I0,DI0,alphaGP,Dmu,Rhomv0,Tmax,1,1,no,sigma,obnoise);
    mv0.track = Trjtemp{1,1};
    mv0.track(:,3) = mv0.track(:,3)*0.1;
    Mtemp;
    clear mv0no Mtemp Trjtemp;
end

[ trjans ] = tracksolverlsp1( moviethis, mv0, sigma, psfdecay, bsize, D, alpha, lambda, deltaI, itemax );

Ithis = mean(trjtrue.track(:,3));

[ vtrue ] = visualizetrjtrack( [trjtrue.track(:,2:-1:1),trjtrue.track(:,3)/Ithis,trjtrue.track(:,4)], S(1), zoom, 1);

[ vans ] = visualizetrjtrack( [trjans.track(:,1:2),trjans.track(:,3)/Ithis,trjans.track(:,4)], S(1), zoom, 1);

lmv = S(1)*zoom;

mvfull = zeros(lmv,lmv*2,Tmax);
for i=1:Tmax
    mvfull(1:lmv,1:lmv,i) = vtrue(:,:,i);
    mvfull(1:lmv,lmv+1:lmv*2,i) = vans(:,:,i);
end

end






