function [w_tilde,phi_tilde,wj_best,jindex, GMModel] = partial_trajectory_mapping(y,yk,plots)
%% Initialisation
T=size(y,1);
K=size(yk,1);

%% Imitation learning
psi=zeros(T,T);
for k=1:T
    for n=1:T
        psi(k,n)=cos(2*pi*(k-1)*(n-1)/T)-1i*sin(2*pi*(k-1)*(n-1)/T);
    end
end
phi_tilde=1/T*(ctranspose(psi));
w_tilde=psi*y;

M=3;
w=[real(w_tilde), imag(w_tilde)];
GMModel = fitgmdist(w,M,'RegularizationValue',0.01);

%% Recomposition test
y_recomposed=phi_tilde*w_tilde;
error=y-real(y_recomposed);
if(plots)
    figure; hold on;
    x=linspace(1,T,T)';
    plot(x,y,'r')
    plot(x,real(y_recomposed),'b')
end

%% Partial trajectory mapping
lambda=1e-8;
phiK_tilde=phi_tilde(1:K,:);
% Calculate the wj best candidates
for j=1:M
    sigma=GMModel.Sigma(:,:,j);
    mu=GMModel.mu(j,:)';

    phiKH_tilde=transpose(real(phiK_tilde))-1i*transpose(imag(phiK_tilde));
    invsigma=inv(sigma);
    invsigma_tilde=invsigma(1,1)+1i*invsigma(2,2);
    mu_tilde=mu(1)+1i*mu(2);

    left=phiKH_tilde*phiK_tilde+lambda*invsigma_tilde;
    right=phiKH_tilde*yk+lambda*invsigma_tilde*mu_tilde;
    wj(:,j)=left\right;
end

%Choose the j best candidate
val=zeros(M,1);
for j=1:M

    weight=GMModel.ComponentProportion(j);
    mu=GMModel.mu(j,:); 
    mu_tilde=mu(1)+1i*mu(2);

    subs=wj(:,j)-mu_tilde;
    val(j)=sum((sqrt(real(yk-phiK_tilde*wj(:,j)).^2+imag(yk-phiK_tilde*wj(:,j)).^2)).^2)-lambda*log(weight)+lambda*sum((sqrt((real(subs)).^2+(imag(subs)).^2)).^2);
    
end
[~,jindex]=min(val);
wj_best=wj(:,jindex);




