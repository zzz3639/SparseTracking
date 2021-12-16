function [ trjanslocations, mvfull, trjansdetail ] = sequentialonoff( moviethis, sig, psfdecay, bsize, D, alpha, lambda, deltaI, sprate, itemax, Kextend, bcmax, trjtrue, live, jobname, mv0 )
%Usage: 
%  inputs:
%    moviethis: movie to be solved, m*n*T
%    sig: Gaussian sigma
%    psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%    bsize: boundary size
%    D: diffusion rate Gaussian sigma
%    alpha: intensity linkage
%    lambda: LSP norm parameter, array of size Tmax+1, lambda(end) is the penalty of terminate
%    deltaI: LSP norm parameter, array of size Tmax
%    sprate: slash value between likelihood and regularized penalty, smaller value gives more sparse result
%    itemax: maximal iteration number
%    Kextend: number of new tracks one previous track generate
%    bcmax: maximal number of frames back extension achieves
%    trjtrue: true answer to be compared with, during the optimization, have no effect solved trajectories
%        use trjtrue=emptytrack() to run the solver without true answer
%    live: 'yes' or 'no', choose 'yes' to show the current result while running.
%    jobname: a string, as tracksolverlspfull writes and reads files from hard disk, 
%             we use jobname to ensure 2 jobs use different files. 
%    mv0: mv0.pic and mv0.no, for first frame, no problem to be ignored.
%  outputs:
%    trjanslocations: molecule locations in all the frames,
%        trjanslocations.no is T*1 array records noise values of all frames
%        trjanslocations.track is a matrix, [x,y,I,T], 2 position values, intensities, frame
%    mvfull: movie, true answer is shown in the left, solved answer is shown in the right
%    trjansdetail: record the detailed information of sequential on off loop

moviethis = moviethis*sprate;
trjtrue.track(:,3) = trjtrue.track(:,3)*sprate;
frame0decay = 1/3; %new emerging molecule gets smaller sparse penalty
Ddecay = 1/3;
zoom = 2;
rng('shuffle');
s1 = size(moviethis,1);
s2 = size(moviethis,2);
S = [s1,s2];
Tmax = size(moviethis,3);
mg = 7;
smol1 = s1-mg-mg;
smol2 = s2-mg-mg;
sarea = smol1*smol2;
Smol = [smol1,smol2,mg,mg];
Intframe0 = 100;
rateintnew = 1/3;
Nframe0 = ceil(sarea*0.5);
if size(lambda,1)==1
    lambda = lambda';
end
if size(deltaI,1)==1
    deltaI = deltaI';
end

if exist('mv0')
else
    mv0no = sum(sum(moviethis(:,:,1)));
    mv0.no = mv0no;
    picframe0 = molacgen(Smol, Nframe0, ones(smol1,smol2), Intframe0, 0, 0);
    mv0.pic = picframe0(:,1:3);
    clear Mtemp Trjtemp trjtemp ;
end

%prepare
N = size(mv0.pic,1);
mv0full.no = mv0.no;
mv0full.track = zeros(N,5);
mv0full.track(:,1:3) = mv0.pic;
mv0full.track(:,5) = [0:N-1]';
%for first frame
lambdasolver = frame0decay*lambda(1)*ones(Nframe0,1);
deltaIsolver = deltaI(1)*ones(Nframe0,1);
[ trjans ] = tracksolverlspfull( moviethis(:,:,1), mv0full, sig, psfdecay, bsize, D, alpha, lambdasolver, deltaIsolver, itemax, jobname );


mvfull = zeros(S(1)*zoom,2*S(2)*zoom,Tmax);
trjanslocations.track = zeros(0,4);
trjanslocations.no = zeros(0,1);
if Tmax<=bcmax
    trjansdetail=cell(1,1);
else
    trjansdetail=cell(Tmax-bcmax,1);
end

% visualization for the true answer(to compare with)
vtrue = visualizetrjtrack( trjtrue.track(:,[1,2,3,4,5]), S, Tmax, zoom, 1);

%sequential on off loop
for t=2:Tmax
    t
    fw = t;
    bc = max(1,t-bcmax);
    md = floor((fw+bc)/2);
    %forward 
        %unify the molecule labels
    trackthis = trjans.track;
    v = find(trackthis(:,4)>=bc-1);
    [trackthis] = molidunique(trackthis(v,:));
    lbnum = max(trackthis(:,5))+1;

        %extend the existing molecules to next frame.
    v = find(trackthis(:,4)==t-2);
    trjend = trackthis(v,:);
    trjend(:,3) = 0.5*trjend(:,3)/Kextend;
    trjendindex = trjend(:,5);
    Nnext = size(trjend,1);
    trjendfull = zeros(0,5);
    for i=1:Nnext
        u = find(trackthis(:,5)==trjendindex(i));
        trackthis(u,3) = 0.5*trackthis(u,3)/Kextend;
        trjendpath = trackthis(u,:);
        trjendpath(:,5) = i-1;
        trjendfull = [trjendfull;trjendpath];
    end
    [vendfull] = find(trjendfull(:,4)==t-2);
    [vend] = find(trackthis(:,4)==t-2);
    
    tracknew = trackthis;
    tracknew = trackmodify(tracknew, vend, 0, D, Ddecay);
    tracknew = trackextend(tracknew, vend, t-1, D); 
    for i=1:Kextend
        tracknewthis = trackmodify(trjendfull, vendfull, lbnum+(i-1)*Nnext, D, Ddecay);
        tracknew = [tracknew;tracknewthis];
    end
     
    for i=1:Kextend-1
        tracknewthis = trackmodify(trjendfull, vendfull, lbnum+(Kextend+i-1)*Nnext, D, Ddecay);
        tracknewthis = trackextend(tracknewthis, vendfull, t-1, D); 
        tracknew = [tracknew;tracknewthis];
    end
        % add possible new molecules to the next frame
    v = find(tracknew(:,4)==t-1);
    if length(v)==0
        Intnovel = Intframe0;
    else
        Intnovel = mean(tracknew(v,3))*rateintnew;
    end
    Rho = ones(smol1,smol2);
    Nnew = ceil(sarea*0.5);
    trackac = molacgen(Smol, Nnew, Rho, Intnovel, t-1, lbnum+(Kextend+Kextend-1)*Nnext);
    tracknew = [tracknew;trackac];

        %save to mv0full and run the forward solver
    mv0full.track = tracknew;
    mv0full.no = [mv0full.no;mv0full.no(end)];
    Nmv0 = size(mv0full.track,1);
    [srange,trange] = molranget(mv0full.track);
    lambdasolver = lambda([trange-srange+1]',1)+lambda(end)*(trange~=t-1);
    deltaIsolver = deltaI([trange-srange+1]',1);
    mv0full.track(:,4) = mv0full.track(:,4) + 1 - bc;
    [ trjans ] = tracksolverlspfull( moviethis(:,:,bc:fw), mv0full, sig, psfdecay, bsize, D, alpha, lambdasolver, deltaIsolver, itemax, jobname );
    
    %backward
        %unify the molecule labels
    trjans.track(:,4) = trjans.track(:,4) + bc -1;
    [trackthis] = molidunique(trjans.track);
    lbnum = max(trackthis(:,5))+1;
    vend = find(trackthis(:,4)==t-1);
    Nlast = sum(vend~=0);
    trjend = trackthis(vend,:);
    trjbcend = trjend;
    trjbcend(:,5) = [0:Nlast-1]';
        %add backward trajectories for each molecule in the last frame
    for i=fw:-1:bc
        trjbc = backinfer(trjbcend, i-1, D, lbnum+(fw-i)*Nlast);
        trackthis = [trackthis;trjbc];
    end
        %save to mv0full and run the backward solver
    mv0full.track = trackthis;
    mv0full.no = trjans.no;
    Nmv0 = size(mv0full.track,1);
    mv0full.track(:,4) = mv0full.track(:,4) - bc +1;
    [srange,trange] = molranget(mv0full.track);
    lambdasolver = lambda([trange-srange+1]',1)+lambda(end)*(trange~=t-1);
    deltaIsolver = deltaI([trange-srange+1]',1);
    [ trjans ] = tracksolverlspfull( moviethis(:,:,bc:fw), mv0full, sig, psfdecay, bsize, D, alpha, lambdasolver, deltaIsolver, itemax, jobname );
    if t>bcmax
        trjansdetail{t-bcmax}=trjans;
    end
    if Tmax<=bcmax
        if t==Tmax
            trjansdetail{1}=trjans;
        end
    end

    %show the last frame this algorithm currently fixed
    Ithis = mean(trjtrue.track(:,3));
    vans = visualizetrjtrack( trjans.track, S, fw-bc+1, zoom, 1);
    vintmax = max(max([vtrue(:,:,t),vans(:,:,end)]))+10;
    mvfull(:,1:S(2)*zoom,1:Tmax) = vtrue(:,:,1:Tmax)/vintmax;
    if Tmax>bcmax
        mdinter = floor((2+bcmax)/2);
        if t==bcmax+1
            mvfull(:,S(2)*zoom+1:2*S(2)*zoom,1:mdinter) = vans(:,:,1:mdinter)/vintmax;
            v = find(trjans.track(:,4)<mdinter);
            trjanslocations.track = [trjanslocations.track;trjans.track(v,1:4)];
            trjanslocations.no = [trjanslocations.no;trjans.no(1:mdinter,1)];
        end
        if t>bcmax+1&&t<Tmax
            mvfull(:,S(2)*zoom+1:2*S(2)*zoom,md) = vans(:,:,mdinter)/vintmax;
            v = find(trjans.track(:,4)==mdinter-1);
            trjnew = trjans.track(v,1:4);
            trjnew(:,4) = trjnew(:,4) + bc - 1;
            trjanslocations.track = [trjanslocations.track;trjnew];
            trjanslocations.no = [trjanslocations.no;trjans.no(mdinter,1)];
        end
        if t==Tmax
            mvfull(:,S(2)*zoom+1:2*S(2)*zoom,md:end) = vans(:,:,mdinter:end)/vintmax;
            v = find(trjans.track(:,4)>=mdinter-1);
            trjnew = trjans.track(v,1:4);
            trjnew(:,4) = trjnew(:,4) + bc - 1;
            trjanslocations.track = [trjanslocations.track;trjnew];
            trjanslocations.no = [trjanslocations.no;trjans.no(mdinter:end,1)];
        end
    else
        if t==Tmax
            mvfull(:,S(2)*zoom+1:2*S(2)*zoom,1:end) = vans(:,:,1:end)/vintmax;
            trjanslocations.track = [trjanslocations.track;trjans.track(:,1:4)];
            trjanslocations.no = [trjanslocations.no;trjans.no];
        end
    end
    if strcmp(live,'yes')
        imshow([vtrue(:,:,md),vans(:,:,floor((1+fw-bc+1)/2))]/vintmax);
        drawnow;
    end
    trjans.track(:,4) = trjans.track(:,4) + bc - 1;
end

%[ trackR, trackP ] = trackingerror1( trjans , trjtrue, IntFilter );
%trackerror = zeros(Tmax,2);
%trackerror(:,1) = (trackR(:,2)+trackP(:,2))./(trackR(:,1)+trackP(:,1));
%trackerror(:,2) = (trackR(:,3)+trackP(:,3))/2;
end

function [trackmdf] = trackmodify(trjendfull, vendfull, idadd, D, Ddecay)
    Nnext = sum(vendfull~=0);
    trackmdf = trjendfull;
    trackmdf(:,5) = trjendfull(:,5)+idadd;
    trackmdf(vendfull,1:2) = trjendfull(vendfull,1:2) + D*Ddecay*randn(Nnext,2);
end

function [trackex] = trackextend(trjendfull,vendfull,t,D)
% extend trajectories survive in frame t to frame t+1
    Nnext = sum(vendfull~=0);
    trjend = trjendfull(vendfull,:);
    trjend(:,1:2) = trjend(:,1:2) + D*randn(Nnext,2);
    trjend(:,4) = t;
    trackex = [trjendfull;trjend];
end

function [trackac] = molacgen(Smol, Nac, Rho, Intnew, t, idadd)
    posac = FNpicgen( Smol, Nac, Rho );
    trackac = [posac,Intnew*ones(Nac,1),t*ones(Nac,1),[0:Nac-1]'+idadd];
end


function [trjbc] = backinfer(trjend, tto, D, idadd)
% back extend particles in frame t to frame tto
    Nlast = size(trjend,1);
    trjbc = trjend;
    if Nlast==0
        return;
    end
    trjbc(:,5) = trjbc(:,5) + idadd;
    trjbc(:,1:2) = trjbc(:,1:2) + D*randn(Nlast,2);
    t = trjend(1,4);
    for j=t-1:-1:tto
        trjbcone = trjbc(1:Nlast,:);
        trjbcone(:,1:2) = trjbcone(1:Nlast,1:2)+D*randn(Nlast,2);
        trjbcone(:,4) = j;
        trjbc = [trjbcone;trjbc];
    end
end


