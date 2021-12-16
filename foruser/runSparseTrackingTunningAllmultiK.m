function [trjans, trjansdetail] = runSparseTrackingTunningAllmultiK(moviethis, movietunning, trjtruetunning, sig, psfdecay, bsize, thint, corenum, jobname, K)
% this function run "runSparseTrackingTunningAll" by multiplying moviethis
%   and movietunning by K. The aim of this function is to avoid sparse penalty and linking
%   penalty dominate the result when the values in moviethis and movietunning are small
% give thint and use trjans, trjansdetail base on the value of moviethis and 
%   movietunning. Multiply and divide are implemented inside this function

%multiply
    trjtruetunning.track(:,3) = trjtruetunning.track(:,3)*K;
    trjtruetunning.no = trjtruetunning.no*K;
    [trjans, trjansdetail] = runSparseTrackingTunningAll(K*moviethis, K*movietunning, trjtruetunning, sig, psfdecay, bsize, K*thint, corenum, jobname);
%divide
    trjans.no = trjans.no/K;
    trjans.track(:,3) = trjans.track(:,3)/K;
    L = length(trjansdetail);
    for i=1:L
        trjtemp = trjansdetail{i};
        trjtemp.no = trjtemp.no/K;
        trjtemp.track(:,3) = trjtemp.track(:,3)/K;
        trjansdetail{i} = trjtemp;
    end
end