function [ moviethis, trjtrue ] = runSimuTracksFunctionPoisson( Smol, N, S, svl, Nframes, Tmax, D, alpha, bgvalue, Imean, Dint )
%RUNSIMUTRACKSFUNCTION Summary of this function goes here
%   Detailed explanation goes here

[Movies, Trj] = SimuTracksPoisson(Smol, S, N/svl, svl-1, Imean, Dint, alpha, D, ones(Smol(1),Smol(2)), zeros(0,3), Nframes, 1, 1, bgvalue*S(1)*S(2), 1, 0);
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

