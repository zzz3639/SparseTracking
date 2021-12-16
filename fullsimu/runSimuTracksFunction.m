function [ moviethis, trjtrue ] = runSimuTracksFunction( Smol, N, S, svl, Nframes, Tmax, D, alpha, bgvalue, Imean, Dint )
% Do the job basically the same to demo 1, simulate molecule motion, activation/disappear and in silicon imaging based on input parameters
% in this function psf_width=1
%   Smol: [s1,s2,margin1,margin2], molecules locate in s1*s2 area, images are (s1+margin1+margin1)*(s2+margin2+margin2)
%   N: average number of molecules in each frame 
%   S: size of images, [s1+margin1+margin1, s2+margin2+margin2]
%   svl: number of frames(in average) a molecule survive(from activated to disappear)
%   Nframes: number of frames simulated
%   Tmax: only last Tmax frames preserved. Expected molecule density could be unbalanced
%         in the frames in the front because of the initialization of frame 1
%   D: motion speed of molecules (gaussian sigma, mean squared displacement=2D^2)
%   alpha: molecule fluorescent intensity fluctuation, log intensity follows gaussian process.
%          function SimuTracks has detailed explanation
%   bgvalue: background fluorescent intensity per pixel
%   Imean, Dint: emitter(molecule) fluorescent intensities follows log normal distribution (LogImean,Dint^2)

[Movies, Trj] = SimuTracks(Smol, S, N/svl, svl-1, Imean, Dint, alpha, D, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, bgvalue*S(1)*S(2), 1, 0);
moviethis = Movies{1,1};
moviethis = moviethis(:,:,end-Tmax+1:end);
trj = Trj{1,1};
v = find(trj(:,4)>=Nframes-Tmax);
trj = trj(v,:);
trj = molidunique(trj);
trj(:,4) = trj(:,4)-Nframes+Tmax;
trjtrue.track = molidunique(trj);
trjtrue.no = bgvalue*S(1)*S(2)*ones(Tmax,1);

end

