clear;
clc;
N = 8;
falconans = cell(N,1);
trjtruefull = cell(N,1);
for i=1:N
    load(['CompareRho',num2str(i)]);
    [molresults, mvresults] = runtrackingfalcon(movietunning);
    falconans{i} = molresults;
    trjtruefull{i} = trjtruetunning;
end
