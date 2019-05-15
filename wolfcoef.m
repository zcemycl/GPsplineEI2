function [ma,mb,Caa,Cbb,Cab] = wolfcoef(T,C,nsr,sig,m)
sizT = size(T);i0 = find(C(:,1)==0);
c1 = 0.05; c2 = 0.8; 
Y = [T(:,2);T(:,3)];
[kTT,kdTT,dkTT,dkdTT] = kernel(T(:,1),T(:,1),sig,m);
Nmat = sig*nsr*diag(ones(1,sizT(1))); 
KTT = [kTT+Nmat,kdTT;...
            dkTT,dkdTT]; 
[kTN,kdTN,dkTN,dkdTN] = kernel(T(:,1),C(:,1),sig,m);

% need
KTTinv = inv(KTT); 
Kvec = [kTN;dkTN]; Kdvec = [kdTN;dkdTN];
knn = zeros(length(C(:,1)),1); 
dknn = zeros(length(C(:,1)),1); 
kdnn = zeros(length(C(:,1)),1); 
dkdnn = zeros(length(C(:,1)),1); 
for i = 1:length(C(:,1))
    [knn(i),kdnn(i),dknn(i),dkdnn(i)]=kernel(C(i,1),...
        C(i,1),sig,m);
end
[~,~,dkn0,~]=kernel(C(:,1),C(i0,1),sig,m);
[k0n,~,dk0n,dkd0n]=kernel(C(i0,1),C(:,1),sig,m);
% gram matrix
GT = Kvec'*KTTinv; GdT = Kdvec'*KTTinv;
% mean mu
mu = GT*Y; dmu = GdT*Y; 
mu0 = mu(i0); dmu0 = dmu(i0);

% ma,mb
ma = mu0-mu+c1*C(:,1)*dmu0;
mb = dmu-c2*dmu0;

% cov k()
cov = knn - diag(GT*Kvec); 
% dcov = dknn - diag(GdT*Kvec);
covd = kdnn - diag(GT*Kdvec);
dcovd=dkdnn - diag(GdT*Kdvec); 
cov00 = knn(i0) - GT(i0,:)*Kvec(:,i0);
covd00 = kdnn(i0) - GT(i0,:)*Kdvec(:,i0);
dcovd00= dkdnn(i0) - GdT(i0,:)*Kdvec(:,i0);
cov0n = [k0n - GT(i0,:)*Kvec]';
dcov0n = [dk0n - GdT(i0,:)*Kvec]';
dcovd0n = [dkd0n - GdT(i0,:)*Kdvec]'; 
dcovn0 = dkn0 - GdT*Kvec(:,i0);

% Caa, Cbb, Cab
Caa=cov00+(c1*C(:,1)).^2*dcovd00+cov+2*(c1*C(:,1).*(covd00-dcov0n)-cov0n);
Cbb=c2^2*covd00-2*c2*dcovd0n+dcovd;
Cab=-c2*(covd00+c1*C(:,1)*dcovd00)+c2*dcov0n+dcovn0+c1*C(:,1).*dcovd0n-covd;


end