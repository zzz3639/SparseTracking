function [jobnumarray]=parajobnum( m, corenum )
% for example, m=8,corenum=3, then jobnumarray=[3;3;2]
jobnumarray = zeros(corenum,1);
jobnumarray(:,1) = floor(m/corenum);
r = mod(m,corenum);
jobnumarray(1:r,1) = jobnumarray(1:r,1) + 1;

end
