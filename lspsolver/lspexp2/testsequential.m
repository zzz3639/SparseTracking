function [ Mvfulls, Trjans ] = testsequential( Movies, Trjtrues, sig, psfdecay, bsize, Dmu, alpha, lambda, deltaI, maxite, Kextend )
%test the ability of sequential solver
%  Usage: [ output_args ] = testsequential( Movies, Trjtrues, sig, bsize, psfdecay, Dmu, alpha, lambda, deltaI, maxite, Kextend )

L = length(Movies);
Mvfulls = cell(L,1);
Trjans = cell(L,1);
for i=1:L
    moviethis = Movies{i};
    trjtrue = Trjtrues{i};
    [ trjans, mvfull ] = sequentialtracking( moviethis, sig, psfdecay, bsize, Dmu, alpha, lambda, deltaI, maxite, Kextend, trjtrue );
    Mvfulls{i,1} = mvfull;
    Trjans{i,1} = trjans;
end

end

