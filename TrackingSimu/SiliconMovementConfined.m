function [Trj] = SiliconMovementConfined(AC,SV,I0,DI0,alpha,Dmu,Centerinfo,S,Rho,Pic0,Tmax,Nrepeat)
%Modified in 2018.7.2 by ZHANG Haowen
%Usage: [Trj]=SiliconMovementCenter(AC,SV,I0,DI0,alpha,Dmu,Pic0,Tmax,Nrepeat)
% this function is basically the same to SiliconMovement, the difference is:
%   there is a circular barrier molecules can not pass, the barrier is defined by Centerinfo
% position of first pixel is (0,0)
% inputs:
%  AC: activation rate, AC*Rho_ij is the number of molecules activated of pixel area i j per frame
%  SV: average survival length of each activated molecule
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alpha: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: diffusion parameter, gaussian sigma/pixels
%  Centerinfo: molecules are confined inside or outside a circle defined by:
%             [centerx,centery,radius]
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
% define output array
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
        Drun = confinedmove(picnSV(:,1:2),Centerinfo,Dmu);
        picnSV(:,1:2) = Drun(:,1:2);
        picnSV(:,3) = exp(alpha*(log(pico(idxSV,3))-LogI0)+sqrt(1-alpha*alpha)*(DI0*randn(length(idxSV),1))+LogI0);
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
        Trjthis = [Trjthis;picnSV;picnAC];
        pico = [picnSV;picnAC];
    end
    %%% engine end
    Trj{l} = Trjthis;
end
end


function [Drun] = confinedmove(locationlist,Centerinfo,Dmu)
    % count inside/outside statues and do Gaussian motion
    radius = Centerinfo(3);
    n = size(locationlist,1);
    D = repmat(Centerinfo(1,1:2),n,1) - locationlist;
    Dl = sqrt(D(:,1).*D(:,1)+D(:,2).*D(:,2));
    Area = (Dl > radius);
    
    locationmove =  locationlist + Dmu*randn(n,2);
    Dmove = repmat(Centerinfo(1,1:2),n,1) - locationmove;
    Dlmove = sqrt(Dmove(:,1).*Dmove(:,1)+Dmove(:,2).*Dmove(:,2));
    Areamove = (Dlmove > radius);   

    Drun = locationmove;
    % outside move to inside
    idx1 = find(Area & (~Areamove));
    % correct outside to inside move
    Dlnew1 = repmat(2*radius,length(idx1),1)- Dlmove(idx1,1);
    Drun(idx1,:) = repmat(Centerinfo(1,1:2),length(idx1),1) - Dmove(idx1,:) .* repmat(Dlnew1./Dlmove(idx1,:),1,2);
    % inside move to outside
    idx2 = find((~Area) & Areamove);
    % correct inside to outside move
    Dlnew2 = repmat(2*radius,length(idx2),1)- Dlmove(idx2,1);
    Drun(idx2,:) = repmat(Centerinfo(1,1:2),length(idx2),1) - Dmove(idx2,:) .* repmat(Dlnew2./Dlmove(idx2,:),1,2);
end

