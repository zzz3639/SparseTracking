function [Trj] = SiliconMovementDecelerate(AC,SV,I0,DI0,alpha,Dmu,S,Rho,Pic0,Tmax,Nrepeat)
%Modified in 2018.7.2 by ZHANG Haowen
%Usage: [Trj]=SiliconMovementDecelerate(AC,SV,I0,DI0,alpha,Dmu,Pic0,Tmax,Nrepeat)
% this function is basically the same to SiliconMovement, the difference is:
%   Dmu here is a matrix instead of a scalar, diffusion rate of all pixel areas are specified 
% position of first pixel is (0,0)
% inputs:
%  AC: activation rate, AC*Rho_ij is the number of molecules activated of pixel area i j per frame
%  SV: average survival length of each activated molecule
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alpha: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: matrix of size s1*s2, diffusion parameter of each pixel area, gaussian sigma/pixels
%  S: [s1,s2], area of focus
%  Rho: The distribution of activated molecules, is normalized inside this function
%       matrix of size s1*s2
%  Pic0: activated molecules in frame 0, [X,Y,I]
%  Tmax: total number of frames simulated, first frame 0, last frame Tmax-1
%  Nrepeat: simulate Nrepeat tracks
% Output:
%  Trj: is a cell array of length Nrepeat, each element is a trajectory.
%      each trajectory is a matrix of 5 columns, [x,y,I,t,id], id is the identification of this molecule, start from 0

%%% Initialization
rng('shuffle');
%define output array
Trj = cell(Nrepeat,1);
% define intermedia values
PSV = 1.0/(SV+1);
LogI0 = log(I0);
% define image positions list
s1 = S(1);
s2 = S(2);
cx1 = [0:s1-1]';
cx1 = repmat(cx1,1,s2);
cx2 = [0:s2-1];
cx2 = repmat(cx2,s1,1);
cx = [reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];
Rho = Rho/sum(sum(Rho));
for l=1:Nrepeat
    %%% engine prepare
    k=size(Pic0,1);
    Trjthis = zeros(size(Pic0,1),5);
    Trjthis(:,1:3) = Pic0;
    Trjthis(:,5) = [0:size(Pic0,1)-1]';
    pico = Trjthis;
    %%% engine start
    for t=1:Tmax-1
        %randomize survived molecules
        idxSV = binornd(1,1-PSV,size(pico,1),1);
        idxSV = find(idxSV);
        %randomize activated molecules
        MAC = poissrnd(AC*Rho);
        NAC = sum(sum(MAC));
        Nthis = length(idxSV) + NAC;
        picnSV = pico(idxSV,:);
        picnSV(:,4) = t;
        picnSV(:,5) = pico(idxSV,5);
        %move the survived molecules
        Dthis = zeros(length(idxSV),1);
        idxSVmatrix = idxpair2matrix(picnSV(:,1:2),S);
        idxSVlist = (idxSVmatrix(:,2)-1)*s1 + idxSVmatrix(:,1);
        Dthis = Dmu(idxSVlist);
        picnSV(:,1) = pico(idxSV,1) + Dthis.*randn(length(idxSV),1);
        picnSV(:,2) = pico(idxSV,2) + Dthis.*randn(length(idxSV),1);
        picnSV(:,3) = exp(alpha*(log(pico(idxSV,3))-LogI0)+sqrt(1-alpha*alpha)*(DI0*randn(length(idxSV),1))+LogI0);
        %generate activated molecules
        picnAC = zeros(NAC,5);
        picnAC(:,4) = t;
        picnAC(:,5) = [k:k+NAC-1]';
        picnAC(:,3) = exp(LogI0+DI0*randn(NAC,1));
        MACline = reshape(MAC,s1*s2,1);
        V = find(MACline);
        j=0;
        for i=1:length(V)
            MACline(V(i));
            picnAC(j+1:j+MACline(V(i)),1:2) = repmat(cx(V(i),:),MACline(V(i)),1) + rand(MACline(V(i)),2) - 0.5;
            j = j+MACline(V(i));
        end
        k = k+NAC;
        %add new molecules in this frame into trajectory list
        Trjthis = [Trjthis;picnSV;picnAC];
        pico = [picnSV;picnAC];
    end
    %%% engine end
    Trj{l} = Trjthis;
end
end

function [idx]=idxpair2matrix(pos,S)
% pos: [x,y], n*2 matrix
% idx: [lineidx,columnidx], n*2 matrix
s1 = S(1);
s2 = S(2);
xgrid = round(pos(:,1));
ygrid = round(pos(:,2));

% for elements in xgrid out of range [0,s2-1], assign 0 if xgrid(i)<0, assign s2-1 if xgrid(i)>s2-1 
xgrid = xgrid - (xgrid<0).*xgrid;
xgrid = xgrid - (xgrid>(s2-1)).*(xgrid-s2+1);

% for elements in ygrid out of range [0,s1-1], assign 0 if ygrid(i)<0, assign s1-1 if ygrid(i)>s1-1 
ygrid = ygrid - (ygrid<0).*ygrid;
ygrid = ygrid - (ygrid>(s1-1)).*(ygrid-s1+1);

% column and line index start from 1
idx = [ygrid+1,xgrid+1];

end



