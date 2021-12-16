function [trjans, trjansdetail] = runSparseTrackingTunningAll(moviethis, movietunning, trjtruetunning, sig, psfdecay, bsize, thint, corenum, jobname)
    [Dfull,Iffull] = trjinfo(trjtruetunning.track,'no');
    alpha = 1/mean(Iffull);
    alpha = min(alpha,20);
    D = sqrt(mean(Dfull.^2)/2);
    [ trjans, trjansdetail ] = runSparseTrackingTunning( moviethis, movietunning, trjtruetunning, sig, psfdecay, bsize, D, alpha, thint, corenum, jobname );
end