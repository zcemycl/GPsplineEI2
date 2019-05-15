% Step samples

% Inputs
% x0 initial vectors, a step size picked by Hessian
% dir direction vector, range rangeafter point 0
% num num of extra samples, noise on objective

% Variable names
% T train: tt training step sample, ot objective train step, 
%          gt gradient of train step (with noise)
% C curve: tc curve step sample, oc objective curve step, 
%          gc gradient curve step (without noise)
% B bar:   two points of the slope bar (from train)
%          (x0,y0),(x1,y1)

function [T,C,B] = stepsamples(dir,a,x0,range,num,noise)
% cost function and its gradient % can be defined externally later
fun = @(x)3*x(1)^2+x(2)^2+55*x(3)^2+2*x(4)^2+x(5)^2;
gra = @(x)2*[3*x(1),x(2),55*x(3),2*x(4),x(5)]';
funa = @(ac)fun(x0+ac*dir);
graa = @(ac)gra(x0+ac*dir)'*dir;

% curve first (since without noise)
tc = [0:range/10000:range]';
oc = fvv(funa, tc);gc = fvv(graa, tc);
C = [tc,oc,gc];

% train sample (with noise)
tt = [0,a]';
te = datasample(tc,num);
tt = cat(1,tt,te);
lt = length(tt); nv = noise*randn(lt,1);
ot = fvv(funa,tt)+nv;gt = fvv(graa,tt);
T = [tt,ot,gt];

% slope bar
B = zeros(lt,4); %[x0,x1,y0,y1]
bl = 0.02*range; % bar length
for i = 1:lt
    ic = ot(i) - gt(i)*tt(i);
    xx = [tt(i)-bl,tt(i)+bl];
    yy = [gt(i)*xx(1)+ic,gt(i)*xx(2)+ic];
    B(i,:) = [xx,yy];
end

end