clear;
clc;
rho = 2;
Nframes = 3000;
Tmax = 2000;
Nframestunning = 300;
Tmaxtunning = 100;
muscale = 7;
Smol = [28,28,7,7];
svl = 20;
D = 0.2;
alpha = 0.8;
bgvalue = 50;
Imean = [500,750,1000,1500,2000];
Dint = 0.5;
Nset = length(Imean);
for i=1:Nset
    imean = Imean(i);
    N = Smol(1)*Smol(2)*rho/muscale/muscale;
    S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
    [ moviethis, trjtrue ] = runSimuTracksFunction( Smol, N, S, svl, Nframes, Tmax, D, alpha, bgvalue, imean, Dint );
    [ movietunning, trjtruetunning ] = runSimuTracksFunction( Smol, N, S, svl, Nframestunning, Tmaxtunning, D, alpha, bgvalue, imean, Dint );
    save(['CompareImean',num2str(imean),'.mat']);
end

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
