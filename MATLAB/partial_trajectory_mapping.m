function [E_mu,w_data,muj,phi,wj_best,jindex, GMModel] = partial_trajectory_mapping(y,yk,plots)
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
A=[eye(T),1i*eye(T)];
phi=phi_tilde*A;
w_tilde=psi*y;

M=10;
w=[real(w_tilde); imag(w_tilde)];
GMModel = fitgmdist(w',M,'RegularizationValue',0.01);

%% Recomposition test
y_recomposed=phi*w;
%error=y-real(y_recomposed);
if(plots)
    figure; hold on;
    x=linspace(1,T,T)';
    plot(x,y(:,1),'r')
    plot(x,real(y_recomposed(:,1)),'b')
end

%% Partial trajectory mapping
lambda=1e-8;
phiK=phi(1:K,:);

% Calculate the wj best candidates
wj=zeros(2*T,M);
for j=1:M
    sigma=GMModel.Sigma(:,:,j);
    mu=GMModel.mu(j,:);
    phiKH = ctranspose(phiK);
    invsigma=eye(size(sigma))/sigma;

    left=phiKH*phiK+lambda*invsigma;
    right=phiKH*yk+lambda*invsigma*mu';
    wj(:,j)=left\right;
end

%Choose the j best candidate
val=zeros(M,1);
for j=1:M

    
    weight=GMModel.ComponentProportion(j);
    mu=GMModel.mu(j,:);
    sigma=GMModel.Sigma(:,:,j);

    subs=wj(:,j)-mu';
    val(j)=sum((sqrt(real(yk-phiK*wj(:,j)).^2+imag(yk-phiK*wj(:,j)).^2)).^2)-lambda*log(weight)+lambda*sum(((sqrt((real(subs)).^2+(imag(subs)).^2)).^2)'/sigma);

end
[~,jindex]=min(val,[],1);
% wj_best=zeros(2*T,1);
% muj=zeros(2*T,1);

wj_best=wj(:,jindex);
muj=GMModel.mu(jindex,:)';

w_data=w;
E_mu=(GMModel.ComponentProportion*GMModel.mu)';

% %plot to test best mu_j*
% figure; hold on;
% plot(real(phi*muj))
% plot(y(:,1))
