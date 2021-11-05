function Y_star = Gauss_Regression(B,Y,X,X_star,param,derivative)
%Solve the equation Y* = Kx*x * inv(Kxx+sigma*I) * Y with(Kxx')=exp(-(x-x')^2/(2*l^2))
%derivative: gives the derivative in fonction of the index number given as
%derivative

if ~exist('derivative','var')
    % third parameter does not exist, so default it to something
    derivative = 0;
end



n=size(B,2);
m=size(X_star,1);
d=size(X_star,2);
% l=sqrt(var(X));
% l=0.2*ones(1,size(X_star,2));
l=param.l;
g=param.g;
if(d==1)
    X_mat_nm=repmat(X,1,m);
    X_star_mat_mn=repmat(X_star,1,n);

    if(not(derivative))
        K_xstar_x=g*exp(-(X_star_mat_mn.'-X_mat_nm).^2./(2.*l.^2));
    else
        K_xstar_x=g*(X_star_mat_mn.'-X_mat_nm).*((-1)./(l.^2)).*exp(-(X_star_mat_mn.'-X_mat_nm).^2./(2.*l.^2));
    end

    Y_star = K_xstar_x.'*B*Y;

else
    X_mat_nm=repmat(X,1,1,m);
    X_star_mat_mn=repmat(X_star,1,1,n);
    K_xstar_x=ones(n,m);

    if(not(derivative))
        for k=1:2
            X_star_mat_mn_1d=reshape(X_star_mat_mn(:,k,:),[m,n]);
            X_mat_nm_1d=reshape(X_mat_nm(:,k,:),[n,m]);
            K_tmp=g*exp(-(X_star_mat_mn_1d(:,:).'-X_mat_nm_1d(:,:)).^2./(2.*l.^2));
            K_xstar_x=K_xstar_x.*K_tmp;
        end
    else
        for k=1:2
            X_star_mat_mn_1d=reshape(X_star_mat_mn(:,k,:),[m,n]);
            X_mat_nm_1d=reshape(X_mat_nm(:,k,:),[n,m]);
            if(k==derivative)
                K_tmp=g*(X_star_mat_mn_1d.'-X_mat_nm_1d).*((-1)./(l.^2)).*exp(-(X_star_mat_mn_1d.'-X_mat_nm_1d).^2./(2.*l.^2));
            else
                K_tmp=g*exp(-(X_star_mat_mn_1d.'-X_mat_nm_1d).^2./(2.*l.^2));
            end
            K_xstar_x=K_xstar_x.*K_tmp;
        end
    end

    Y_star = K_xstar_x'*B*Y;
end