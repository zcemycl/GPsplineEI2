clear all; clc;
% test script
% settings
% num: number of extra steps
% noise: noise level
% sig: signal level
% nsr: noise-to-signal ratio
num = 1; noise = 0.5;
nsr = 20; sig = 20;
m = 20;

% input 
pt4 = [0.8561,4.6657,-0.0115,-2.0839,3.6095]';
dir = [0.2064,0.1797,0.0114,0.4667,0.1391]'; 
x0 = pt4; a = 0.3840;

% samples
[T,C,B] = stepsamples(dir,a,x0,5*a,num,noise);

% posterior/predictive distribution
% (negative log marginal likelihood
% for hyparameter estimation built-in)
[mu,cov,nsr,sig,newC,pt] = posterior(T,C,nsr,sig,m); 

% expected improvement
u = EI(mu,cov); [~,indm] = max(u);
fprintf('expected improvement: %f \n',C(indm,1))

% wolfe termination probability
% [newC,pt] = probWolfe(T,C,nsr,sig,m);

% plot stuff
subplot(3,1,1)
grid on
plotmucov(mu,cov,C(:,1)); hold on;
plotsamples(T,C,B);
legend('confidence','\mu','original curve',...
    'train data','Location','bestoutside')

subplot(3,1,2)
plot(C(:,1),u)
hold on; grid on
xlim([0,max(C(:,1))])
legend('improvement','Location',...
    'bestoutside');
ylabel('expected improvement')

subplot(3,1,3)
plot(newC,pt)
hold on; grid on
xlim([0,max(C(:,1))])
legend('termination','Location',...
    'bestoutside');
ylabel('Wolfe probability')
