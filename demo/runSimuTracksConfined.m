rho = 2; %molecule density, mol/mum^2
muscale = 7; % pixels/mum
Smol = [28,28,7,7];
N = Smol(1)*Smol(2)*rho/muscale/muscale;
S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
svl = 20;
zoom = 10;
Nframes = 6000; % Nframes frames are simulated, only the last Tmax frames are preserved
Tmax = 5000;
Dmu = 0.2;
center = [Smol(1)/2-0.5,Smol(2)/2-0.5];
Centerinfo = [center,7];
[Movies, Trj] = SimuTracksConfined(Smol, S, N/svl, svl-1, 1000, 0.5, 0.7, Dmu, Centerinfo, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, 50*S(1)*S(2), 1, 0);
moviethis = Movies{1,1};
moviethis = moviethis(:,:,end-Tmax+1:end);
trj = Trj{1,1};
v = find(trj(:,4)>=Nframes-Tmax);
trj = trj(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframes+Tmax;
trjtrue.track = molidunique(trj);
trjtrue.no = 50*S(1)*S(2)*ones(Tmax,1);
%mvtrue = visualizetrjtrack( trjtrue.track, S, 10, 1 );
%scale bar for 2*pixel(2*psfsigma)
%mvtrue(3:4 ,zoom+1:zoom+2*zoom,:) = max(max(max(mvtrue)));

% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
lambda = 0.5*[20:19+Tmax,0];
deltaI = 10*ones(1,Tmax);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
