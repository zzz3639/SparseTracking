function [ trjans, mvfull ] = sequentialtracking( moviethis, sig, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, Kextend, trjtrue, mv0 )
%Usage: 
%    moviethis: movie to be solved, m*n*T
%    sig: Gaussian sigma
%    psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%    bsize: boundary size
%    D: diffusion rate Gaussian sigma
%    alpha: intensity linkage
%    lambda: LSP norm parameter
%    deltaI: LSP norm parameter
%    itemax: maximal iteration number
%    Kextend: number of new tracks one previous track generate
%    mv0: mv0.pic and mv0.no, for first frame
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
if exist('mv0')
else
    mv0no = sum(sum(moviethis(:,:,1)));
    mv0.no = mv0no;
    Rhomv0 = ones(smol1,smol2)/(smol1*smol2);
    Rhomv0 = ceil(sarea/2)*Rhomv0;
    [Mtemp, Trjtemp] = SimuTracksNaive(Smol,S,1000,0.5,0.7,0.4,Rhomv0,1,1,1,0,sig,0);
    trjtemp = Trjtemp{1,1};
    mv0.pic = [trjtemp(:,2:-1:1),trjtemp(:,3)/10];
    clear Mtemp Trjtemp trjtemp ;
end

N = size(mv0.pic,1);
mv0full.no = mv0.no;
mv0full.track = zeros(N,5);
mv0full.track(:,1:3) = mv0.pic;
mv0full.track(:,5) = [0:N-1]';
for t=1:Tmax
    %run EM algorithm
    [ trjans ] = tracksolverlsp1( moviethis(:,:,1:t), mv0full, sig, psfdecay, bsize, D, alpha, lambda, deltaI*t, itemax );
    Ithis = mean(trjtrue.track(:,3));
    vtrue = visualizetrjtrack( trjtrue.track(:,[2,1,3,4,5]), 40, 10, 1);
    vans = visualizetrjtrack( trjans.track, 40, 10, 1);
    if t==Tmax
        mvfull = cat(2,vtrue,vans)/Ithis;
    end
    vintmax = max(max([vtrue(:,:,t),vans(:,:,end)]));
    imshow([vtrue(:,:,t),vans(:,:,end)]/vintmax);
    drawnow;
    %unify the molecule labels
    trackthis = trjans.track;
    lbindex = unique(trackthis(:,5));
    lbnum = length(lbindex);
    k = 0;
    for i=1:lbnum
        u = find(trjans.track(:,5)==lbindex(i));
        trackthis(u,5) = i-1;
    end
    %extend the existing molecules to next frame.
    v = find(trackthis(:,4)==t-1);
    trjend = trackthis(v,:);
    trjend(:,3) = trjend(:,3)/Kextend;
    trjendindex = trjend(:,5);
    Nnext = size(trjend,1);
    trjendfull = zeros(0,5);
    for i=1:Nnext
        u = find(trackthis(:,5)==trjendindex(i));
        trackthis(u,3) = trackthis(u,3)/Kextend;
        trjendfull = [trjendfull;trackthis(u,:)];
    end
    for i=1:Kextend-1
        trackthis = [trackthis;[trjendfull(:,1:4),trjendfull(:,5)+lbnum+(i-1)*Nnext]];
    end
    
    trackframenext = trjend;
    trackframenext(:,1:2) = trjend(:,1:2) + D*randn(Nnext,2);
    trackframenext(:,4) = trjend(:,4) + 1;
    trackthis = [trackthis;trackframenext];  
    for i=1:Kextend-1
        trackframenext = trjend;
        trackframenext(:,4) = trackframenext(:,4) + 1;
        trackframenext(:,5) = trackframenext(:,5) + lbnum + (i-1)*Nnext;
        trackframenext(:,1:2) = trjend(:,1:2) + D*randn(Nnext,2);
        trackthis = [trackthis;trackframenext];
    end
    trackthis(:,1:2) = trackthis(:,2:-1:1);
    mv0full.track = trackthis;
    mv0full.no = [mv0full.no;mv0full.no(end)];
end

end


