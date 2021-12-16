function [ Movies ] = SiliconImagingTracks( m, s, trj, no, sigma, obnoise )
%Modified in 2016.11.29 by ZHANG Haowen
% On silicon imaging of some moving particles.
% Usage: [ Movies ] = SiliconTracks( m, s, trj, no, sigma, obnoise )
%    m: number of repeats
%    s=[s1,s2], field of vision, all start from 0
%    trj=[X,Y,Intensity,T,id], molecule list, all start from 0
%    no: noise, count by photon number
%    sigma: PSF width
%    obnoise: Additional noise, for each pixel, signal=photoncount+sqrt(photoncount)*randn()*obnoise

%%% initialize
rng('shuffle');
Tmax = max(trj(:,4))+1;
s1=s(1);
s2=s(2);
[U,V] = sort(trj(:,4));
trj = trj(V,:);
NT = zeros(Tmax+1,1);
Sp = sparse(U+1,[1:length(U)]',ones(length(U),1));
Sp = sum(Sp,2);
NT(2:end,:) = full(Sp);
NT = cumsum(NT);
% define image positions list
cx1=[0:s1-1]';
cx1=repmat(cx1,1,s2);
cx2=[0:s2-1];
cx2=repmat(cx2,s1,1);
cx=[reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];
% define output variable
Movies=cell(m,1);

%%% generate poisson images
for i=1:m
    moviethis = zeros(s1,s2,Tmax);
    for t=1:Tmax
        pic = trj(NT(t)+1:NT(t+1),1:3);
        n = NT(t+1) - NT(t);
        img=zeros(s1*s2,1);
        for j=1:n
            PSFdis = PSF(pic(j,1:2),cx,sigma);
            imgdif = [mnrnd(round(pic(j,3)),[PSFdis;abs(1-sum(PSFdis))])]';
            img=img+imgdif(1:end-1,:);
        end
        img=img+[Noise(s1,s2,round(no))]';
        img=reshape(img,s1,s2);
        moviethis(:,:,t) = img;
    end
    Movies{i} = moviethis;
end

%%% add additional noise
for i=1:m
    moviethis = Movies{i};
    for t=1:Tmax
        img = moviethis(:,:,t);
        imgb = zeros(size(img));
        imgb = img+obnoise*randn(s1,s2).*sqrt(img);
        imgb = imgb-(imgb<0).*imgb;
        moviethis(:,:,t) = imgb;
    end
    Movies{i} = moviethis;
end

end

function b2=Noise(s1,s2,PhotonNum)
    b2=mnrnd(PhotonNum,repmat(1/s1/s2,s1*s2,1));
end

%PSF function
function V=PSF(mu,cx,sig)
    V=normpdf(cx(:,1),mu(1),sig).*normpdf(cx(:,2),mu(2),sig);
    S = sum(sum(V));
    if S>1
        V = V/S;
    end
end


