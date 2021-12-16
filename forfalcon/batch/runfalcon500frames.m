clear;
clc;
P = '/home/zzz/trackingexperiments/datasetbac/';
N = 8;
thint = 200;
falconans = cell(N,1);
trjtruefull = cell(N,1);
for i=1:N
    load([P,'CompareRho',num2str(i),'.mat']);
    trjfalcon = runfalcontunning(  moviethis, movietunning, trjtruetunning, thint );
    falconans{i} = trjfalcon;
    trjtruefull{i} = trjtrue;
end

save('falconres.mat');