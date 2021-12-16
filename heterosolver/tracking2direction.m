function [ trjanslocations, mvfull, trjansdetail, tempmvfullforeback, temptrjdetailforeback ] = tracking2direction( moviethis, sig, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, Kextend, bcmax, trjtrue, jobname )
%TRACKING2DIRECTION run sequentialonoff 1toT and Tto1, then combine the results of the 2 directions into a better solution
%   Detailed explanation goes here

if size(lambda,1)==1
    lambda = lambda';
end
if size(deltaI,1)==1
    deltaI = deltaI';
end
trjtrueback = trjtrue;
trjtrueback.no = trjtrue.no(end:-1:1);
T = size(moviethis,3);
zoom = 10;
s1 = size(moviethis,1);
s2 = size(moviethis,2);
S = [s1,s2];
trjtrueback.track(:,4) = T-1-trjtrue.track(:,4);
trjtrueforeback = cell(2,1);
trjtrueforeback{1} = trjtrue;
trjtrueforeback{2} = trjtrueback;

movieforeback = cell(2,1);
movieforeback{1} = moviethis;
movieforeback{2} = moviethis(:,:,end:-1:1);

%run sequentialonoff forward and backward
temptrjdetailforeback=cell(2,1);
tempmvfullforeback=cell(2,1);
for i=1:2
    [ temptrjanslocations, tempmvfull, temptrjdetail ] = sequentialonoff ... 
        (movieforeback{i}, sig, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, Kextend, bcmax, trjtrueforeback{i}, 'yes', jobname);
    temptrjdetailforeback{i} = temptrjdetail;
    tempmvfullforeback{i} = tempmvfull;
end

%combine the results and output final results
L = length(temptrjdetailforeback{1});
trjansdetail = cell(L,1);

if bcmax>=T
    %L==1 in this case
    trjthisfore = temptrjdetailforeback{1}{1};
    trjthisback = temptrjdetailforeback{2}{1};
    trjthisback.no = trjthisback.no(end:-1:1);
    trjthisback.track(:,4) = T - 1 - trjthisback.track(:,4);
    trjthisback.track = sorttrack(trjthisback.track);
    trjunion = trackunion(trjthisfore,trjthisback);
else
    for i=1:L
        i
        %frame i to frame i+bcmax
        trjthisfore = temptrjdetailforeback{1}{i};
        trjthisback = temptrjdetailforeback{2}{L-i+1};
        moviefragment = moviethis(:,:,i:i+bcmax);
        %reverse the backward trajectories
        trjthisback.no = trjthisback.no(end:-1:1);
        trjthisback.track(:,4) = bcmax - trjthisback.track(:,4);
        trjthisback.track = sorttrack(trjthisback.track);
        %merge forward trajectory set with backward trajectory set
        trjunion = trackunion(trjthisfore,trjthisback);
        %run the solver
        [srange,trange] = molranget(trjunion.track);
        lambdasolver = lambda([trange-srange+1]',1);
        deltaIsolver = deltaI([trange-srange+1]',1);
        [ trjans ] = tracksolverlspfull( moviefragment, trjunion, sig, psfdecay, bsize, D, alpha, lambdasolver, deltaIsolver, itemax, jobname );
        trjansdetail{i} = trjans;
    end
    % output trjanslocations and mvfull
    trjanslocations.track = zeros(0,4);
    trjanslocations.no = zeros(0,1);
    mvfull = zeros(S(1)*zoom,2*S(2)*zoom,T);
    mdidx = floor((2+bcmax)/2);
    vtrue = visualizetrjtrack( trjtrue.track, S, zoom, 1 );
    mvfull(:,1:S(2)*zoom,1:T) = vtrue(:,:,1:T);
    trjans = trjansdetail{1};
    vans = visualizetrjtrack( trjans.track, S, zoom, 1 );
    mvfull(:,S(2)*zoom+1:2*S(2)*zoom,1:mdidx-1) = vans(:,:,1:mdidx-1);
    for i=1:mdidx-1
        v = find(trjans.track(:,4)==i-1);
        trjanslocations.track = [trjanslocations.track;trjans.track(v,1:4)];
        trjanslocations.no = [trjanslocations.no;trjans.no(i)];
    end
    for i=1:L
        trjans = trjansdetail{i};
        vans = visualizetrjtrack( trjans.track, S, zoom, 1 );
        mvfull(:,S(2)*zoom+1:2*S(2)*zoom,i-1+mdidx) = vans(:,:,mdidx);
        v = find(trjans.track(:,4)==mdidx-1);
        trjlocationsthis = trjans.track(v,1:4);
        trjlocationsthis(:,4) = i-2+mdidx;
        trjanslocations.track = [trjanslocations.track;trjlocationsthis];
        trjanslocations.no = [trjanslocations.no;trjans.no(mdidx)];
    end
    trjans = trjansdetail{end};
    vans = visualizetrjtrack( trjans.track, S, zoom, 1 );
    mdbcnum = bcmax - mdidx;
    mvfull(:,S(2)*zoom+1:2*S(2)*zoom,end-mdbcnum:end) = vans(:,:,mdidx+1:end);
    for i=mdidx+1:bcmax+1
        v = find(trjans.track(:,4)==i-1);
        trjlocationsthis = trjans.track(v,1:4);
        trjlocationsthis(:,4) = i+T-bcmax-2;
        trjanslocations.track = [trjanslocations.track;trjlocationsthis];
        trjanslocations.no = [trjanslocations.no;trjans.no(i)];
    end
end

end
