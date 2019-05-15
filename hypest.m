function [newnsr,newsig,nlml] = hypest(T,Y,nsr,sig,m)
leni1 = length(T);
[KA,KO] = Kf(T,sig,m);

Nmat = nsr*diag(ones(1,leni1)); 
F = [Nmat,zeros(leni1,leni1);...
    zeros(leni1,leni1),zeros(leni1,leni1)];
Ky = KO+F; Zyi = KA+KO+F;
Zy = inv(Zyi); Z = Zy/sig; z = Z*Y;
B = Z-z*z'; 
dKyd = 2*nsr*F; 
ddKyd = 4*nsr*F;
sp = Y'*Zy*Y/(2*leni1-3); 
dsp = -Y'*Zy*dKyd*Zy*Y;
dKc = sp*dKyd+dsp*Ky;
ddsp = 2*Y'*Zy*dKyd*Zy*dKyd*Zy*Y/(2*leni1-3)-Y'*Zy*ddKyd*Y;

term1 = 0.5*trace(dKc*Z*dKc*(2*z*z'-Z));
term2 = 2*nsr*dsp*trace(B*F);
term3 = 2*nsr*sp*trace(B*F);
term4 = 0.5*ddsp*trace(B*Ky);

% marginal likelihood
nlml = 0.5*Y'*Z*Y+0.5*log(det(sig*Zyi))+log(2*pi)*leni1/2;

% gradient and Hessian
g = 0.5*trace(B*dKc);
H = term1+term2+term3+term4; 

% new parameters
newnsr = nsr - exp(2*g/H);

[KAn,KOn] = Kf(T,sig,m);
Nmatn = newnsr*diag(ones(1,leni1)); 
Fn = [Nmatn,zeros(leni1,leni1);...
    zeros(leni1,leni1),zeros(leni1,leni1)];
Zyin = KAn+KOn+Fn; Zyn = inv(Zyin);
newsig = Y'*Zyn*Y/(2*leni1-3);

end