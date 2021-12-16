function [ trjans ] = tracksolverlspfull( moviethis, mv0, sig, psfdecay, bsize, D, alpha, lambda, deltaI, itemax, jobname )
%  Usage: [ trjans ] = tracksolverlsp1( moviethis, mv0, sig, psfdecay, bsize, D, alpha, lambda, deltaI, itemax )
%  call afulllsp.out, execute EM solver
%    moviethis: movie to be solved, m*n*T
%    mv0: mv0.track and mv0.no
%    sig: Gaussian sigma
%    psfdecay: range(pixels) of psf, out of psfdecay psf are set to 0
%    bsize: boundary size
%    D: diffusion rate
%    alpha: intensity linkage
%    lambda: LSP norm parameter, array of length, size(mv0.track,1)
%    deltaI: LSP norm parameter, array of length, size(mv0.track,1)
%    itemax: maximal iteration number
f1 = 'paratemp';
f2 = 'movietemp';
f3 = 'mv0temp';
f4 = 'anstemp';
f1 = [jobname,f1];
f2 = [jobname,f2];
f3 = [jobname,f3];
f4 = [jobname,f4];
% if the maximal intensity of a trajectory is less than molfiltervalue, remove this trajectory
molfiltervalue = 50.0;

fp1 = fopen(f1, 'w');
fprintf(fp1,'%d\n%d\n',[bsize,psfdecay]);
fprintf(fp1,'%f\n',sig);
fprintf(fp1,'%f\n%f\n',[D,alpha]);
fprintf(fp1,'%d\n',itemax);
fprintf(fp1,'%f\n',molfiltervalue);
fclose(fp1);
if size(lambda,1)==1
    dlmwrite(f1,lambda,'-append','delimiter',' ');
else
    dlmwrite(f1,lambda','-append','delimiter',' ');
end
if size(deltaI,1)==1
    dlmwrite(f1,deltaI,'-append','delimiter',' ');
else
    dlmwrite(f1,deltaI','-append','delimiter',' ');
end

fp2 = fopen(f2, 'w');
fprintf(fp2,'%d %d %d\n',[size(moviethis,1),size(moviethis,2),size(moviethis,3)]);
fclose(fp2);
for i=1:size(moviethis,3)
    dlmwrite(f2,moviethis(:,:,i),'-append','delimiter',' ');
end

fp3 = fopen(f3, 'w');
fprintf(fp3,'%d %d\n',[size(mv0.track,1),size(moviethis,3)]);
fclose(fp3);
dlmwrite(f3,mv0.track,'-append','delimiter',' ');
if size(mv0.no,1)==1
    dlmwrite(f3,[mv0.no],'-append','delimiter',' ');
else
    dlmwrite(f3,[mv0.no]','-append','delimiter',' ');
end

system(['afulllsp.out ', f2, ' ', f3, ' <', f1, ' >', f4]);
A = importdata(f4);
ns = A.data(end,1);
trjans.track = A.data(1:end-1-size(moviethis,3),:);
trjans.track(:,1:2) = trjans.track(:,2:-1:1);
trjans.no = A.data(end-size(moviethis,3):end-1,1);

end





