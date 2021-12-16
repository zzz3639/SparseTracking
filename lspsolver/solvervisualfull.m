function [trjans, mvfull, mv0] = solvervisualfull(moviethis, trjtrue, sigma, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, mv0)
%Usage: [trjans, mvfull] = solvervisual(moviethis, trjtrue, sigma, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, mv0)
%  moviethis: movie to be solved, m*n*T
%  trjtrue: the true answer
%  sigma: PSF width
%  psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%  bsize: boundary size
%  D: diffusion rate
%  alpha: intensity linkage
%  lambda: LSP norm parameter, array of length, size(mv0.track,1)
%  deltaI: LSP norm parameter, array of length, size(mv0.track,1)
%  itemax: maximal iteration number
%  mv0: mv0.track and mv0.no

zoom = 20;
mg = 7;
s1 = size(moviethis,1);
s2 = size(moviethis,2);
Tmax = size(moviethis,3);
smol1 = s1-mg-mg;
smol2 = s2-mg-mg;
sarea = smol1*smol2;
Smol = [smol1,smol2,mg,mg];
if exist('mv0')
else
    mv0no = sum(sum(moviethis));
    mv0.no = zeros(Tmax,1);
    mv0.no(:,1) = mv0no(1,1,:);
    Rhomv0 = ones(smol1,smol2)/(smol1*smol2);
    Rhomv0 = ceil(sarea/2)*Rhomv0;
    [Mtemp, Trjtemp] = SimuTracksNaive(Smol,[s1,s2],1000,0.5,0.7,0.4,Rhomv0,Tmax,1,1,50,sigma,0);
    mv0.track = Trjtemp{1,1};
    mv0.track(:,3) = mv0.track(:,3)*0.1;
    Mtemp;
    clear mv0no Mtemp Trjtemp;
end

mv0.track(:,1:2) = mv0.track(:,2:-1:1);
Nmol = max(mv0.track(:,5)) + 1;
Mmol = sparse(mv0.track(:,4)+1,mv0.track(:,5)+1,ones(size(mv0.track,1),1));
Lmol = full(sum(Mmol,1));

if length(lambda)==1
    lambda = lambda*ones(Nmol,1);
end

if length(deltaI)==1
    deltaI = deltaI*Lmol';
end

[ trjans ] = tracksolverlspfull( moviethis, mv0, sigma, psfdecay, bsize, D, alpha, lambda, deltaI, itemax );

Ithis = mean(trjtrue.track(:,3));

[ vtrue ] = visualizetrjtrack( [trjtrue.track(:,2:-1:1),trjtrue.track(:,3)/Ithis,trjtrue.track(:,4)], s1, zoom, 1);

[ vans ] = visualizetrjtrack( [trjans.track(:,1:2),trjans.track(:,3)/Ithis,trjans.track(:,4)], s1, zoom, 1);

lmv = s1*zoom;

mvfull = zeros(lmv,lmv*2,Tmax);
for i=1:Tmax
    mvfull(1:lmv,1:lmv,i) = vtrue(:,:,i);
    mvfull(1:lmv,lmv+1:lmv*2,i) = vans(:,:,i);
end

end






