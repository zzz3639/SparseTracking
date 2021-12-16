function [ alignfrom, alignto ] = emitteralias( emitterfrom, emitterto, D, alpha, th )
% align emitter list "emitterfrom" to emitter list "emitterto"
% emitterlists are matrix [x,y,I]
% distance is defined as ((x1-x2)^2+(y1-y2)^2)/D^2 + abs(log(I1)-log(I2))/alpha
% only distance<th is valid align
% alignfrom and alignto are line numbers, where emitter alignfrom(i) align to alignto(i)

m1 = size(emitterfrom,1);
m2 = size(emitterto,1);
distmatrix = zeros(m1,m2);
for i=1:m2
    distmatrix(:,i) = ((emitterfrom(:,1)-emitterto(i,1)).^2+(emitterfrom(:,2)-emitterto(i,2)).^2)/D^2 + ...
        abs(log(emitterfrom(:,3))-log(emitterto(i,3)))/alpha;
end

[u,v] = min(distmatrix,[],2);

valid = find(u<th);
valid = valid;

alignfrom = valid;
alignto = v(valid);

end


