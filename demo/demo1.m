
% emitter(molecule) density, average number every 7*7 pixels area, approximate to number of molecules per um^2
rho = 2;
% length measurement, 7 pixel=1um
muscale = 7;
% molecules locate in 28*28 area, images are (28+7+7)*(28+7+7) as there's 7 pixel margin in each direction
Smol = [28,28,7,7];
% average number of molecule in each frame 
N = Smol(1)*Smol(2)*rho/muscale/muscale;
% size of images
S = [Smol(1)+2*Smol(3),Smol(2)+2*Smol(4)];
% number of frames(in average) a molecule survive(from activated to disappear)
svl = 10;
% number of frames simulated
Nframes = 1000;
% only last Tmax frames preserved. Expected molecule density could be unbalanced
% in the frames in the front because of the initialization of frame 1
Tmax = 100;
% motion speed of molecules (gaussian sigma, mean squared displacement=2D^2)
D = 0.2;
% background fluorescent intensity per pixel
bgpixel = 50;
% psf width (Gaussian sigma)
sig = 1;
% emitter(molecule) fluorescent intensities follows log normal distribution (Log(Imean),Dint^2)
Imean = 1000;
Dint = 0.5;
% molecule fluorescent intensity fluctuation, log intensity simulated by gaussian process.
%     in function SimuTracks there's detailed explanation
alpha = 0.8;
% do the simulation, simulate molecule motion, activation/disappear and in silicon imaging 
[Movies, Trj] = SimuTracks(Smol, S, N/svl, svl-1, Imean, Dint, alpha, D, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, bgpixel*S(1)*S(2), sig, 0);
% extract the answer(s)
moviethis = Movies{1,1};
trj = Trj{1,1};
% Keep the last Tmax frames, remove the others
moviethis = moviethis(:,:,end-Tmax+1:end);
v = find(trj(:,4)>=Nframes-Tmax);
trj = trj(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframes+Tmax;
trjtrue.track = molidunique(trj);
trjtrue.no = 50*S(1)*S(1)*ones(Tmax,1);

%Visualization of simulated images and ground truth
% for visualization, images and ground truth are both 210*210 in size in visualization,210=42*5
zoom = 5;
% vtrj = visualizetrjtrack( trj(:,1:4), S(1), zoom , 1);
% implay(vtrj);
%[ trjans, mvfull, trjansdetail ] = sequentialonoff( moviethis(:,:,:), 1, 5, 0, D, 3.3, lambda, deltaI, 5000, 2, 10, trjtrue, 'yes', '' );
