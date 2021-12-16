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
%old
%GPalpha = [0.6,0.7,0.8,0.9,0.95]; 
%gaussian process alpha, log(It) = gpalpha*log(It-1) + sqrt(1-gpalpha^2)*randn(...)
%new
GPalpha = [0.4,0.5]; 
bgvalue = 50;
Imean = 1000;
Dint = 0.5;
Nset = length(GPalpha);
for i=1:Nset
    gpalpha = GPalpha(i);
    N = Smol(1)*Smol(2)*rho/muscale/muscale;
    S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
    [ moviethis, trjtrue ] = runSimuTracksFunction( Smol, N, S, svl, Nframes, Tmax, D, gpalpha, bgvalue, Imean, Dint );
    [ movietunning, trjtruetunning ] = runSimuTracksFunction( Smol, N, S, svl, Nframestunning, Tmaxtunning, D, gpalpha, bgvalue, Imean, Dint );
    save(['CompareGPalpha',num2str(gpalpha),'.mat']);
end

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
