function u = EI(mu,cov)
n = min(mu);
term1 = (n-mu).*(1+erf((n-mu)./sqrt(2*cov)))/2;
term2 = sqrt(cov/2/pi).*exp(-(n-mu)/2./cov);

u = term1 + term2;

end