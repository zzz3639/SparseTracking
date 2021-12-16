function [Fchecked,Tchecked] = checkmapping(F, T)
%original mapping is F to T, remove contradict mappings(two or more in F map to one in T).
%F, array of positive integers, elements are unique
%T, array of positive integers, same size as F, may not be unique

l = length(F);
sp = sparse(F,T,ones(l,1));
w = full(sum(sp,1));
v = find(w(T)==1);
Fchecked = [F(v)];
if size(Fchecked,1)==1
    Fchecked = Fchecked';
end
Tchecked = [T(v)];
if size(Tchecked,1)==1
    Tchecked = Tchecked';
end

end

