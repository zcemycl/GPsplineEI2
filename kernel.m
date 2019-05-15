function [k,kd,dk,dkd] = kernel(t,t_,sig,m)
leni1 = length(t); leni2 = length(t_);
k = zeros(leni1,leni2);kd = zeros(leni1,leni2);
dk = zeros(leni1,leni2);dkd = zeros(leni1,leni2);

M = [m^3/3,m^2/2;m^2/2,m];

for i = 1:leni1
    x = t(i);
    h = [1,0;x,1]; 
    for j = 1:leni2
        y = t_(j);
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
        Ka = h'*M*h_;
        k(i,j) = sig*(ko+Ka(1,1));
        kd(i,j) = sig*(kdo+Ka(1,2));
        dk(i,j) = sig*(dko+Ka(2,1));
        dkd(i,j) = sig*(dkdo+Ka(2,2));
    end
    
end

end