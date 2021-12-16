function [Trj] = SiliconMovementNaive(I0,DI0,alpha,Dmu,S,Pic0,Tmax,Nrepeat)
%Modified in 2017.11.06 by ZHANG Haowen
%molecule movement without activation or bleeching
%Usage: [Trj] = SiliconMovementNaive(I0,DI0,alpha,Dmu,S,Pic0,Tmax,Nrepeat)
%  I0,DI0: Parameters of Intensity distribution (log normal distribution(LogI0,DI0^2))
%  alpha: the intensity of a certain molecule over time is considered to follow a gaussian process: I(m,t+1)=alpha*I(m,t)+sqrt(1-alpha*alpha)*randn()
%  Dmu: diffusion parameter, normal, sigma/pixels
%  S: [s1,s2], area of focus
%  Pic0: activated molecules in frame 0, [X,Y,I]
%  Tmax: total number of frames simulated, first frame 0, last frame Tmax-1
%  Nrepeat: simulate Nrepeat tracks
%  Output is a cell array of length Nrepeat, each element is a trajectory.
%  each trajectory is a matrix of 5 columns, [x,y,I,t,id], id is the identification of this molecule, start from 0

%%% Initialization
rng('shuffle');
%define output array
Trj = cell(Nrepeat,1);
% define intermedia values
LogI0 = log(I0);
% define image positions list
s1 = S(1);
s2 = S(2);
cx1 = [0:s1-1]';
cx1 = repmat(cx1,1,s2);
cx2 = [0:s2-1];
cx2 = repmat(cx2,s1,1);
cx = [reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];
for l=1:Nrepeat
    %%% engine prepare
    k=size(Pic0,1);
    Trjthis = zeros(size(Pic0,1),5);
    Trjthis(:,1:3) = Pic0;
    Trjthis(:,5) = [0:size(Pic0,1)-1]';
    pico = Trjthis;
    %%% engine start
    for t=1:Tmax-1
        idxSV = ones(size(pico,1),1);
        idxSV = find(idxSV);
        picnSV = pico(idxSV,:);
        picnSV(:,4) = t;
        picnSV(:,5) = pico(idxSV,5);
        picnSV(:,1) = pico(idxSV,1) + Dmu*randn(length(idxSV),1);
        picnSV(:,2) = pico(idxSV,2) + Dmu*randn(length(idxSV),1);
        picnSV(:,3) = exp(alpha*(log(pico(idxSV,3))-LogI0)+sqrt(1-alpha*alpha)*(DI0*randn(length(idxSV),1))+LogI0);
        Trjthis = [Trjthis;picnSV];
        pico = [picnSV];
    end
    %%% engine end
    Trj{l} = Trjthis;
end
end


