function [trjcut] = trjtruncate(trj,ts,tt)
% trajectories from ts to tt are extracted

v = find( trj(:,4)<=tt & trj(:,4)>=ts );
trjcut = trj(v,:);

end
