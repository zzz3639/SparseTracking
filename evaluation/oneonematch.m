function [ match1, match2 ] = oneonematch( trj1, trj2 )
%
%   trj1,trj2: [x,y,I,t]

Tmin = 0;
Tmax = max(max(trj1(:,4)), max(trj2(:,4)));
match1 = zeros(size(trj1,1),1);
match2 = zeros(size(trj2,1),1);

for t=Tmin:Tmax
    v1 = find(trj1(:,4)==t);
    v2 = find(trj2(:,4)==t);
    n1 = length(v1);
    n2 = length(v2);
    if n1==0||n2==0
        continue;
    end
    % particles in frame t;
    trjthis1 = trj1(v1,1:3);
    trjthis2 = trj2(v2,1:3);
    % particles aligned to nearest neighbor 
    align1 = MolListAlign(trjthis1, trjthis2);
    align2 = MolListAlign(trjthis2, trjthis1);
    % find one on one match
    M1 = zeros(n1,n2);
    M2 = zeros(n1,n2);
    for i=1:n1
        M1(i,align1(i,7)) = 1;
    end
    for i=1:n2
        M2(align2(i,7),i) = 1;
    end
    for i=1:n1
        j = align1(i,7);
        if M2(i,j)==1&&sum(M1(:,j))==1&&sum(M2(i,:))==1
            % valid one one match
            match1(v1(i)) = v2(j);
            match2(v2(j)) = v1(i);
        end
    end
    
end


end

