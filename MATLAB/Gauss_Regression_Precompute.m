function B = Gauss_Regression_Precompute(X,param)
%Return the constant part of the function Y* = Kx*x * inv(Kxx+sigma*I) * Y with(Kxx')=exp(-(x-x')^2/(2*l^2))
%So B = inv(Kxx+sigma*I) * Y


n=size(X,1);
d=size(X,2);

sigma=param.sigma;
l=param.l;
g=param.g;

if(d==1)
    X_mat_nn=repmat(X,[ones(1,size(X,2)),n]);
    K_xx=g*exp(-(X_mat_nn.'-X_mat_nn).^2./(2.*l.^2));
    identity = eye(size(K_xx));
    K_xx_noisy=K_xx+sigma*identity;
    B=inv(K_xx_noisy);

else
    X_mat_nn=repmat(X,1,1,n);
    K_xx=ones(size(X,1));
    for k=1:2
        X_mat_1d=reshape(X_mat_nn(:,k,:),[n,n]);
        k_tmp=g*exp(-(X_mat_1d(:,:).'-X_mat_1d(:,:)).^2./(2.*l.^2));
        K_xx=K_xx.*k_tmp;

    end
    identity = eye(size(K_xx));
    K_xx_noisy=K_xx+sigma*identity;
    B = inv(K_xx_noisy);
end