function [Dave, Ifave] = trjansdetailinfo( trjansdetail, thinit, noprint )
%TRJANSDETAILINFO Summary of this function goes here
%   Detailed explanation goes here
N = length(trjansdetail);
Dfull = cell(N,1);
Iffull = cell(N,1);
for i=1:N
    trjansthis = trjansdetail{i};
    trjfilted = intensityfilter5column( trjansthis.track, thinit );
    [Dthis, Ifthis] = trjinfo(trjfilted,'no');
    Dfull{i} = Dthis;
    Iffull{i} = Ifthis;
end

Dfullarray = zeros(0,1);
Iffullarray = zeros(0,1);
for i=1:N
    Dfullarray = [Dfullarray;Dfull{i}];
    Iffullarray = [Iffullarray;Iffull{i}];
end

Dave = sqrt(mean(Dfullarray.^2));
Ifave = mean(Iffullarray);

if ~exist('noprint')
    fprintf('\n Average squared displacement: %f', Dave);
    fprintf('\n median displacement: %f', median(Dfullarray));
    fprintf('\n Intensity fluctuation: %f', Ifave);
    fprintf('\n');
end

end

