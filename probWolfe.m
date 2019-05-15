function [newsamplets,pt] = probWolfe(T,C,nsr,sig,m)
[ma,mb,Caa,Cbb,Cab] = wolfcoef(T,C,nsr,sig,m);
[alim,blim,rhot] = coefcdf(ma,mb,Caa,Cbb,Cab);
% alim=real(alim);
% blim=real(blim);
% rhot=real(rhot);
% 
% disp(isreal(alim))
% disp(isreal(blim))
% disp(isreal(rhot))
% disp(size(alim))
% disp(size(blim))
% disp(size(rhot))

% just check the real no.
Aind=find(alim==real(alim));
Bind=find(blim==real(blim));
Rind=find(rhot==real(rhot));
ind = union(union(Aind,Bind),Rind);
% disp(blim(length(ind)))
% disp(ind(length(ind)))
blim = blim(ind); alim = alim(ind); rhot = rhot(ind);
samplets = C(:,1);
newsamplets = samplets(ind);
pt = zeros(1,length(blim));
for i = 1:length(blim)
    pt(i) = bvn(alim(i),inf,blim(i),inf,rhot(i));
end

end
