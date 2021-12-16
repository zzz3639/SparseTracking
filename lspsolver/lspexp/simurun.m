clear all
Smol = [26,26,7,7];
S = [40,40];
maxite = 5000;
Ntry = 5;
Tmax = 10;
Dmu = 0.3;
sig = 1;
bsize = 0;
psfdecay = 5;
alpha = 5.0;
lambda = 20;
deltaI = 1000;
movies = cell(Ntry,1);
trjtrues = cell(Ntry,1);
trjanses = cell(Ntry,2);
mvfulls = cell(Ntry,2);
for i=1:Ntry
    i
    [moviethis, trjtrue, trjans, mvfull] = simusolvervisual(Smol,S,1000,0.5,0.7,Dmu,28*ones(26,26)/26/26,Tmax,S(1)*S(2)*50,sig,0, psfdecay, bsize, Dmu, alpha, lambda, deltaI, maxite);
    [trjans2, mvfull2] = solvervisual(moviethis, trjtrue, sig, psfdecay, bsize, Dmu, alpha, lambda, deltaI, maxite);
    movies{i,1} = moviethis;
    trjtrues{i,1} = trjtrue;
    trjanses{i,1} = trjans;
    trjanses{i,2} = trjans2;
    mvfulls{i,1} = mvfull;
    mvfulls{i,2} = mvfull2;
end