function [ Result ] = MolListAlign(Pic, PicTrue)
% Usage: [ Result ] = MolListAlign(Pic, PicTrue)
%   Pic is the estimated molecule list, in which first 2 columns record positions
%   PicTrue is the true answer

if size(PicTrue,1)==0
    Result = zeros(size(Pic,1),7);
    Result(:,1:3) = Pic;
    Result(:,4:6) = 1000;
    Result(:,7) = 0;
    return;
end
D=dist(Pic(:,1:2),PicTrue(:,1:2));
[u,v]=min(D,[],1);
Result=[Pic(:,1:2),Pic(:,3),(min(D,[],1))',Pic(:,1:2)-PicTrue(v,1:2),v'];

end

function D=dist(mu1,mu2)
    m=size(mu1,1);
    n=size(mu2,1);
    D=zeros(n,m);
    for i=1:m
        T=mu2-repmat(mu1(i,:),n,1);
        D(:,i)=T(:,1).*T(:,1)+T(:,2).*T(:,2);
    end
    D=sqrt(D);
end

