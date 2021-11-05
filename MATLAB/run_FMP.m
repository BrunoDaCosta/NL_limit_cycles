%% Setting Data
close all
clear

plots=true;
load('datatest12.mat', 'traj');
Kini=20;
Kmax=size(traj.data,1);
Xk=traj.data(Kini:Kmax,1);
Xy=traj.data(Kini:Kmax,2);
K=Kmax-Kini+1;


% w_tild = fft(Xk);
% Xk_ = ifft(w_tild);
% w=[real(w_tild),imag(w_tild)];
% k=2;
% GMModel = fitgmdist(w,k,'RegularizationValue',0.01);

psi=zeros(K,K);
for k=0:K-1
    for n=0:K-1
        psi(k+1,n+1)=cos((2*pi*k*n)/K)-i*sin((2*pi*k*n)/K);
    end
end

w_tilde=psi'*Xk;

phi_tilde=(1/K)*ctranspose(psi);
phi=zeros(K,K,2);
for k=1:K
    for n=1:K
        phi(k,n,:)=[real(phi_tilde(n,k)),imag(phi_tilde(n,k))];
    end
end


N=2;
w=[real(w_tilde), imag(w_tilde)];
GMModel = fitgmdist(w,N,'RegularizationValue',0.01);



lambda=1e-8;
for n=1:N
    %phi_=(phi)

    %psiH=inv(ctranspose(phi);
    %val(n)=(psiH()*psi+lambda*inv(GMModel.Sigma(:,:,n))))*(ctranspose(psi)*Xk+lambda*inv(GMModel.Sigma(:,:,n))+GMModel.mu(n,:))
end

%t=ctranspose(psi)(N,:)*psi+lambda*inv(GMModel.Sigma(:,:,n));

if(plots)
    figure;
    plot(Xk);

    figure;
    y = [zeros(1000,1);ones(1000,1)];
    h = gscatter(w(:,1),w(:,2));
    hold on
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
    g = gca;
    fcontour(gmPDF,[g.XLim g.YLim])
    title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
    legend(h,'datapoints')
    hold off
end

