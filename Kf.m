function [KA,KO] = Kf(T,sig,m)
leni1 = length(T);
Ka = zeros(leni1,leni1); 
Kad = zeros(leni1,leni1);
dKa = zeros(leni1,leni1);
dKad = zeros(leni1,leni1);
Ko = zeros(leni1,leni1); 
Kod = zeros(leni1,leni1);
dKo = zeros(leni1,leni1);
dKod = zeros(leni1,leni1);

M = [m^3/3,m^2/2;m^2/2,m];

for i = 1:leni1
    x = T(i);
    h = [1,0;x,1]; 
    for j = 1:leni1
        y = T(j);
        h_ = [1,0;y,1]; 
        if x>y
            ko = (y)^3/3+0.5*(x-y)*(y)^2;
            kdo = x*y-0.5*y^2;
            dko = 0.5*y^2;
            dkdo = y;
        elseif x<y
            ko = (x)^3/3+0.5*(y-x)*(x)^2;
            kdo = 0.5*x^2;
            dko = x*y-0.5*x^2;
            dkdo = x;
        elseif x==y
            ko = (x)^3/3;
            kdo = x*y-0.5*y^2;
            dko = x*y-0.5*x^2;
            dkdo = x;
        end  
        ka = h'*M*h_;
        Ka(i,j) = ka(1,1);
        Kad(i,j) = ka(1,2);
        dKa(i,j) = ka(2,1);
        dKad(i,j) = ka(2,2);
        Ko(i,j) = ko;
        Kod(i,j) = kdo;
        dKo(i,j) = dko;
        dKod(i,j) = dkdo;
    end  
end

KA = [Ka,Kad;dKa,dKad];
KO = [Ko,Kod;dKo,dKod];

end