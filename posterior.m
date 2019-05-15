% Input
% T train set samples
% C curve set samples
function [mu,cov,nsr,sig,newsamplets,pt] = posterior(T,C,nsr,sig,m)
sizT = size(T); 

% marginal likelihood for hyperparameters
% nsr and signal intensity
NLML = []; j = 0; NSR = []; SIG = [];
while true
    j = j+1; 
    NSR = [NSR,nsr]; SIG = [SIG,sig];
    [nsr,sig,nlml] = hypest(T(:,1),[T(:,2);T(:,3)],nsr,sig,m); 
    fprintf('iteration %d nlml: %f \n',[j,nlml])
    NLML = [NLML,nlml];
    if j > 1
        if ~isreal(nlml)
            break
        end
        if NLML(j)-NLML(j-1)>.01
            break
        end
        if j>10000
            break
        end
    end 
end
[~,ind] = min(NLML); ad = 0;
while true
    sig = SIG(ind-ad); nsr = NSR(ind-ad);

    % train data matrix 2Tx2T
    [kTT,kdTT,dkTT,dkdTT] = kernel(T(:,1),T(:,1),sig,m);
    Nmat = sig*nsr*diag(ones(1,sizT(1))); 
    KTT = [kTT+Nmat,kdTT;...
                dkTT,dkdTT]; 
    KTTinv = inv(KTT); 

    % kvec with curve and train data 2T x N
    [kTN,~,dkTN,~] = kernel(T(:,1),C(:,1),sig,m);
    Kvec = [kTN;dkTN]; 

    % observation vector
    Y = [T(:,2);T(:,3)]; GT = Kvec'*KTTinv;
    
    knn = zeros(length(C(:,1)),1); 
    for i = 1:length(C(:,1))
        knn(i) = kernel(C(i,1),C(i,1),sig,m);
    end

    mu = GT*Y;
    cov = knn - diag(GT*Kvec); 
    
    if isreal(sqrt(cov))
        fprintf('iteration %d hyp is selected!\n',ind-ad)
%         disp(size(GT)); disp(size(Kvec))
        break
    end
    ad = ad + 1;
end
while true
    sig = SIG(ind-ad); nsr = NSR(ind-ad);
    try
        [newsamplets,pt] = probWolfe(T,C,nsr,sig,m);
        fprintf('iteration %d hyp is selected!\n',ind-ad)
        break
    catch
        ad = ad + 1;
    end
end


end