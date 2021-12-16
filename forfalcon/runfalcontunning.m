function [ molresults ] = runfalcontunning(  moviethis, movietunning, trjtruetunning, thint )
% Usage: [ molresults ] = runfalcontunning(  moviethis, movietunning, trjtruetunning, thint )
%   Detailed explanation goes here

[ spbest ] = falconparaselect( movietunning, trjtruetunning, thint );
[molresults] = runtrackingfalcon(moviethis,spbest);

end

