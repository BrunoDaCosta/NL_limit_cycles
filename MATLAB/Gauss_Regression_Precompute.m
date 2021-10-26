function B = Gauss_Regression_Precompute(X)
%Return the constant part of the function Y* = Kx*x * inv(Kxx+sigma*I) * Y with(Kxx')=exp(-(x-x')^2/(2*l^2))
%So B = inv(Kxx+sigma*I) * Y


n=size(X,1);
d=size(X,2);
l=1;
sigma=0.001;

if(d>1)
    v_perm=[3,2,1];
else
    v_perm=[2,1];
end

X_mat_nn=repmat(X,[ones(1,size(X,2)),n]);

K_xx=exp(-(permute(X_mat_nn,v_perm)-X_mat_nn).^2./(2.*l.^2));

if(d>1)
    v_perm=[1,3,2];
else
    v_perm=[1,2];
end
identity=permute(repmat(eye(size(K_xx,1)),[ones(1,min(d,2)),min(d,2)]),v_perm);

K_xx_noisy=K_xx+sigma*identity;
invK_xx=squeeze(zeros(n,min(2,d),n));
if(d>1)
    for k=1:min(2,d)
        invK_xx(:,k,:)=inv(squeeze(K_xx_noisy(:,k,:)));
    end
else

    invK_xx(:,:)=inv(squeeze(K_xx_noisy(:,:)));
end

B = reshape(invK_xx,n*min(d,2),[]);
