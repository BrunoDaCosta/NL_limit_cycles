function Y_star = Gauss_Regression(B,Y,X,X_star,derivative)
%Solve the equation Y* = Kx*x * inv(Kxx+sigma*I) * Y with(Kxx')=exp(-(x-x')^2/(2*l^2))

if ~exist('derivative','var')
    % third parameter does not exist, so default it to something
    derivative = false;
end



n=size(B,2);
m=size(X_star,1);
d=size(X_star,2);
l=1;

if(d>1)
    v_perm=[3,2,1];
else
    v_perm=[2,1];
end

X_mat_nm=repmat(X,[ones(1,size(X,2)),m]);
X_star_mat_mn=repmat(X_star,[ones(1,min(size(X_star,2),2)),n]);


if(not(derivative))
    K_xstar_x=exp(-(permute(X_star_mat_mn,v_perm)-X_mat_nm).^2./(2.*l.^2));
else
    K_xstar_x=((-1)./(l.^2))*exp(-(permute(X_star_mat_mn,v_perm)-X_mat_nm).^2./(2.*l.^2)).*((permute(X_star_mat_mn,v_perm)-X_mat_nm));
end

Y_star = reshape(K_xstar_x,n*min(d,2),[])'*B*Y;
