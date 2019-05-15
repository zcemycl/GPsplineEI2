% plot mean,cov confidence
function plotmucov(mu,cov,samples)
f = [mu-2*sqrt(cov); ...
    flip(mu+2*sqrt(cov),1)];
hold on; fill([samples; flip(samples)], ...
    f, [7.5 7.5 7.5]/8);
plot(samples,mu,'k-');
end