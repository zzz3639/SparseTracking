function [] = directedmotion( ndata, nmol, dratio )
%D Summary of this function goes here
%   Detailed explanation goes here
rng('shuffle');
arraydata=zeros(ndata,2);
for i=1:ndata
    A=repmat([dratio,0],nmol,1)+randn(nmol,2);
    arraydata(i,:)=sum(A);
end

plot(arraydata(:,1),arraydata(:,2),'.');

end

