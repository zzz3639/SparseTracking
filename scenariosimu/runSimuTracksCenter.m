clear;
clc;
rho = 2; %molecule density, mol/mum^2
muscale = 7; % pixels/mum
Smol = [42,42,7,7];
N = Smol(1)*Smol(2)*rho/muscale/muscale;
S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
svl = 20;
Imean = 1000;
Dint = 0.5;
Nframes = 6000; % Nframes frames are simulated, only the last Tmax frames are preserved
Tmax = 5000;
Nframestunning = 300;
Tmaxtunning = 100;
Dmu = 0.3;
alpha = 0.8;
no = 50;
center = [Smol(1)/2-0.5,Smol(2)/2-0.5];
Centerinfo = [center,0.05]; % center location, orbitting movement speed
% generate data with ground truth
[Movies, Trj] = SimuTracksCenter(Smol, S, N/svl, svl-1, Imean, Dint, alpha, Dmu, Centerinfo, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, no*S(1)*S(2), 1, 0);
moviethis = Movies{1,1};
moviethis = moviethis(:,:,end-Tmax+1:end);
trj = Trj{1,1};
v = find(trj(:,4)>=Nframes-Tmax);
trj = trj(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframes+Tmax;
trjtrue.track = molidunique(trj);
trjtrue.no = no*S(1)*S(2)*ones(Tmax,1);
% generate data for hyper-parameter tunning
[Movies, Trj] = SimuTracksCenter(Smol, S, N/svl, svl-1, Imean, Dint, alpha, Dmu, Centerinfo, ones(Smol(1),Smol(2)), zeros(0,3), Nframestunning, 1, 1, no*S(1)*S(2), 1, 0);
movietunning = Movies{1,1};
movietunning = movietunning(:,:,end-Tmaxtunning+1:end);
trj = Trj{1,1};
v = find(trj(:,4)>=Nframestunning-Tmaxtunning);
trj = trj(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframestunning+Tmaxtunning;
trjtruetunning.track = molidunique(trj);
trjtruetunning.no = no*S(1)*S(2)*ones(Tmaxtunning,1);

clear Movies Trj trj;
save(['ScenarioOrbitD03F5000.mat']);

%mvtrue = visualizetrjtrack( trjtrue.track, S, 10, 1 );
%scale bar for 2*pixel(2*psfsigma)
%mvtrue(3:4 ,zoom+1:zoom+2*zoom,:) = max(max(max(mvtrue)));

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
