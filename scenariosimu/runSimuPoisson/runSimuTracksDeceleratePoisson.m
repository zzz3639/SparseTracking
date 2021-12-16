clear;
clc;
rho = 2; %molecule density, mol/mum^2
muscale = 7; % pixels/mum
Smol = [42,42,7,7];
N = Smol(1)*Smol(2)*rho/muscale/muscale;
S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
svl = 20;
Imean = 300;
Dint = 0.3;
alpha = 1;
noscale = 20;
sig = 1;
% Nframes number of trajectories simulated, last Tmax frames preserved
Nframes = 2000;
Tmax = 1000;
Nframestunning = 300;
Tmaxtunning = 100;
Dmufast = 0.3;
Dmuslow = sqrt(Dmufast^2/2);
Dmu = zeros(Smol(1),Smol(2));
center = [Smol(1)/2-0.5,Smol(2)/2-0.5];
radius = 2*muscale;
Dmap = zeros(42,42);
for i=0:Smol(1)-1
    for j=0:Smol(2)-1
        if ((i-center(1))^2 + (j-center(2))^2) < radius^2
            Dmu(i+1,j+1) = Dmuslow;
        else
            Dmu(i+1,j+1) = Dmufast;
            Dmap(i+1,j+1) = 1;
        end
    end
end
% generate trajectory data with correct answer
[Movies, Trj] = SimuTracksDeceleratePoisson(Smol, S, N/svl, svl-1, Imean, Dint, alpha, Dmu, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, noscale*S(1)*S(2), sig, 0);
moviethis = Movies{1,1};
moviethis = moviethis(:,:,end-Tmax+1:end);
trjlong = Trj{1,1};
v = find(trjlong(:,4)>=Nframes-Tmax);
trj = trjlong(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframes+Tmax;
trjtrue.track = molidunique(trj);
trjtrue.no = noscale*S(1)*S(2)*ones(Tmax,1);

% generate trajectory data for hyper-parameter tunning
[Moviestunning, Trjtunning] = SimuTracksDeceleratePoisson(Smol, S, N/svl, svl-1, Imean, Dint, alpha, Dmu, ones(Smol(1),Smol(2)), zeros(0,3), Nframestunning, 1, 1, noscale*S(1)*S(2), sig, 0);
movietunning = Moviestunning{1,1};
movietunning = movietunning(:,:,end-Tmaxtunning+1:end);
trjlong = Trjtunning{1,1};
v = find(trjlong(:,4)>=Nframestunning-Tmaxtunning);
trj = trjlong(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframestunning+Tmaxtunning;
trjtruetunning.track = molidunique(trj);
trjtruetunning.no = noscale*S(1)*S(2)*ones(Tmax,1);
clear trj trjlong Moviestunning Trjtunning Movies Trj
save(['ScenarioDecelerateD03Poisson.mat']);
%mvtrue = visualizetrjtrack( trjtrue.track, S, 10, 1 );
% scale bar for 2*pixel(2*psfsigma)
%mvtrue(3:4 ,zoom+1:zoom+2*zoom,:) = max(max(max(mvtrue)));

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
