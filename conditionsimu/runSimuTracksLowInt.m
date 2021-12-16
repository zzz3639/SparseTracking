clear;
clc;
Rho = [0.5,1,1.5,2,2.5,3,3.5,4];
muscale = 7;
Smol = [28,28,7,7];
S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
svl = 20;
Nframes = 3000;
Tmax = 2000;
Nframestunning = 300;
Tmaxtunning = 100;
D = 0.3;
alpha = 1;
bgvalue = 20;
Imean = 300;
Dint = 0.3;
for i=1:length(Rho)
    rho = Rho(i);
    N = Smol(1)*Smol(2)*rho/muscale/muscale;
    [ moviethis, trjtrue ] = runSimuTracksFunctionPoisson( Smol, N, S, svl, Nframes, Tmax, D, alpha, bgvalue, Imean, Dint );
    [ movietunning, trjtruetunning ] = runSimuTracksFunctionPoisson( Smol, N, S, svl, Nframestunning, Tmaxtunning, D, alpha, bgvalue, Imean, Dint );
    save(['TestLowIntRho',num2str(i),'.mat']);
end

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
