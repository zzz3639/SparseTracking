 function [ pic,no,mv, Time ] = EMcrossnframesboundarysparse( bm,n,sigma,iteration, t , bsize, bdecay, lambda, mv0)
%CROSS2FRAMES Summary of this function goes here
%   This is the primary implementation of sparse linking EM algorithm, slow and doesn't support on/off switch of emitters.
%initialize
tic;
m=size(bm,2);
rng();
s=sqrt(size(bm,1));
%initiate at random
%pic=[rand(n,2)*(s-1)+bsize,ones(n,1)/(n+1)];
%pic=[repmat(pic(:,1:2),1,m),pic(:,3)];
%no=1/(n+1);

%initiate from input
pic=mv0.pic;
no=mv0.no;

cx1=[0:s-1+2*bsize]';
cx=repmat(cx1,1,s+2*bsize);
cx=[reshape(cx',(s+2*bsize)*(s+2*bsize),1),reshape(cx,(s+2*bsize)*(s+2*bsize),1)];

cx1=[bsize:s-1+bsize]';
cxo=repmat(cx1,1,s);
cxo=[reshape(cxo',s*s,1),reshape(cxo,s*s,1)];

sbs=(s+2*bsize)*(s+2*bsize);

W=cell(m,1);
for i=1:m
    W{i} = sparse(zeros(sbs,n));
end

Wn=zeros(sbs,m);

A =sparse(zeros(sbs,n));
Ao=sparse(zeros(s*s,n));

SA=zeros(sbs,m);

SW=zeros(n,m);

S=zeros(1,1);

O=zeros(n,m);

T=zeros(n,m);

D=zeros(m,m);

DD=sparse(zeros(n,n));

P=zeros(m,n);

Zero=10^-10;

bbm=zeros(sbs,m);
for i=1:m
    P1=zeros(s+2*bsize,s+2*bsize);
    P2=reshape(bm(:,i),s,s);
    P1(bsize+1:bsize+s,bsize+1:bsize+s)=P2;
    bbm(:,i)=reshape(P1,(s+2*bsize)*(s+2*bsize),1);
end
bm=bbm;

mv.pic=zeros(n,2*m+1,iteration);
mv.no=zeros(iteration,1);

for i=1:iteration
%E-step
    for mi=1:m
        A=PSF(0,s-1+2*bsize,0,s-1+2*bsize,pic(:,2*mi-1:2*mi),sigma,bdecay);
        Ao=PSF(bsize,s-1+bsize,bsize,s-1+bsize,pic(:,2*mi-1:2*mi),sigma,bdecay);
        SA(:,mi)=A*pic(:,2*m+1)+Noise(cx)*no;
        
        W{mi}= spdiags(ones(sbs,1)./SA(:,mi),0,sbs,sbs) * A * spdiags(pic(:,2*m+1),0,n,n);
        Wn(:,mi)=Noise(cx)*no./SA(:,mi);
        Judge=(cx(:,1)>=bsize).*(cx(:,1)<bsize+s).*(cx(:,2)>=bsize).*(cx(:,2)<bsize+s);
        scale=sum(bm(:,mi))/(sum(Ao,1)*pic(:,2*m+1)+no*s*s/(s+2*bsize)/(s+2*bsize));
        bbm(:,mi)=scale*SA(:,mi).*(1-Judge)+bm(:,mi).*Judge;
    end
%m-step

    for mi=1:m
        SW(:,mi)=(W{mi})'*bbm(:,mi);
    end
    np=sum(sum(bbm));
    no0=no;
    no=sum(sum(Wn.*bbm))+np*lambda*no0;
    pic(:,2*m+1)=sum(SW,2);
    S=no+sum(pic(:,2*m+1));
    no=no/S;
    pic(:,2*m+1)=pic(:,2*m+1)/S;
    
    %first axes
    for mi=1:m
        T(:,mi)=(W{mi})'*sparse(bbm(:,mi).*cx(:,1));
        O(:,mi)=(W{mi})'*sparse(bbm(:,mi));
    end
    D=zeros(m,m);
    for j=1:n
        for mi=1:m
            D(mi,mi)=O(j,mi)+2*t;
        end
        for mi=1:m-1
            D(mi,mi+1)=-t;
            D(mi+1,mi)=-t;
        end
        D(1,1)=D(1,1)-t;
        D(m,m)=D(m,m)-t;
        if rcond(D)>Zero
            P(:,j)=D^(-1)*(T(j,:))';
        end
    end
    for mi=1:m
        pic(:,2*mi-1)=(P(mi,:))';
    end
    
    %second axes
    for mi=1:m
        T(:,mi)=(W{mi})'*sparse(bbm(:,mi).*cx(:,2));
        O(:,mi)=(W{mi})'*sparse(bbm(:,mi));
    end
    D=zeros(m,m);
    for j=1:n
        for mi=1:m
            D(mi,mi)=O(j,mi)+2*t;
        end
        for mi=1:m-1
            D(mi,mi+1)=-t;
            D(mi+1,mi)=-t;
        end
        D(1,1)=D(1,1)-t;
        D(m,m)=D(m,m)-t;
        if rcond(D)>Zero
            P(:,j)=D^(-1)*(T(j,:))';
        end
    end
    for mi=1:m
        pic(:,2*mi)=(P(mi,:))';
    end
        
    mv.no(i)=no;
    mv.pic(:,:,i)=pic;

end

Time=toc

end


function Bp=PSF(L1, U1, L2, U2, mu,sigma,bdecay)
    n1=U1-L1+1;
    n2=U2-L2+1;
    nspots=size(mu,1);
    %Bp=sparse(n1*n2,1);
    I=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    J=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    V=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    num=0;
    for k=1:nspots
        m1=floor(mu(k,1));
        m2=floor(mu(k,2));
        l1=max(m1-bdecay,L1);
        u1=min(m1+bdecay+1,U1);
        l2=max(m2-bdecay,L2);
        u2=min(m2+bdecay+1,U2);
        t1=u1-l1+1;
        t2=u2-l2+1;
        if t1<0
            t1=0;
        end
        if t2<0
            t2=0;
        end
    %for i=l1:u1
    %    Bp((i-L1)*n2+l2-L2+1:(i-L1)*n2+u2-L2+1,1)=1/(2*pi*sigma*sigma)*exp((-(i-mu(1))^2-([l2:u2]'-mu(2)).^2)/(2*sigma*sigma));
    %end
        Bp1=zeros(1,u1-l1+1);
        Bp1(1,1:u1-l1+1)=exp(-(([l1:u1]-mu(k,1)).^2)/(2*sigma*sigma));
        Bp2=zeros(u2-l2+1,1);
        Bp2(1:u2-l2+1,1)=exp(-(([l2:u2]'-mu(k,2)).^2)/(2*sigma*sigma));
        V(num+1:num+t1*t2)=reshape([(Bp2*Bp1)/(2*pi*sigma*sigma)]',t1*t2,1);
        II=(repmat([l1:u1]',1,t2)-L1)*n2+repmat([l2:u2],t1,1)-L2+1;
        I(num+1:num+t1*t2)=reshape(II,t1*t2,1);
        J(num+1:num+t1*t2)=k;
        num=num+t1*t2;
    %ones(n1,1)*ones(1,n2);
    end
    Bp=sparse(I(1:num),J(1:num),V(1:num),n1*n2,nspots);
end

function Bn=Noise(X)
    Bn=ones(size(X,1),1)/size(X,1);
end





