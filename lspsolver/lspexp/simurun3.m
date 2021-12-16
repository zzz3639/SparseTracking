clear all
Smol = [26,26,7,7];
S = [40,40];
maxite = 5000;
Ntry = 5;
movies = cell(Ntry,1);
trjtrues = cell(Ntry,1);
trjanses = cell(Ntry,2);
mvfulls = cell(Ntry,2);
for i=1:Ntry
    i
    [moviethis, trjtrue, trjans, mvfull] = simusolvervisual(Smol,S,1000,0.5,0.7,0.3,42*ones(26,26)/26/26,10,40*40*50,1,0, 5, 0, 0.3, 5.0, 20, 1000, maxite);
    [trjans2, mvfull2] = solvervisual(moviethis, trjtrue, 1, 5, 0, 0.3, 5.0, 20, 1000, maxite);
    movies{i,1} = moviethis;
    trjtrues{i,1} = trjtrue;
    trjanses{i,1} = trjans;
    trjanses{i,2} = trjans2;
    mvfulls{i,1} = mvfull;
    mvfulls{i,2} = mvfull2;
end