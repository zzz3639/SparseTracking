% simulate low density fluorescent image series to test MTT
clear;
clc;
Rho = [0.1,0.2,0.3,0.4,0.5];
for i=1:length(Rho)
    rho = Rho(i);
    muscale = 7;
    Smol = [28,28,7,7];
    N = Smol(1)*Smol(2)*rho/muscale/muscale;
    S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
    svl = 20;
    Nframes = 10000;
    Tmax = 5000;
    Nframestunning = 200;
    Tmaxtunning = 100;
    D = 0.2;
    alpha = 0.8;
    bgvalue = 50;
    Imean = 1000;
    Dint = 0.5;
    [ moviethis, trjtrue ] = runSimuTracksFunction( Smol, N, S, svl, Nframes, Tmax, D, alpha, bgvalue, Imean, Dint );
    [ movietunning, trjtruetunning ] = runSimuTracksFunction( Smol, N, S, svl, Nframestunning, Tmaxtunning, D, alpha, bgvalue, Imean, Dint );
    save(['CompareRhoMTT',num2str(i),'.mat']);
end

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
