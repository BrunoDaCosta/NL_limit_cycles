clear; close all;

%% Creation of dataset
plots=false;

T=1000;
K=120;
x=linspace(0,20*pi,T)';
y=sin(x);

% Same initial values or another one
%yk=sin(x(1:K));
yk=cos(2*x(1:K));

%% partial trajectory mapping 
[w_tilde,phi_tilde,wt,jindex,GMModel]=partial_trajectory_mapping(y,yk,plots);
ynew=zeros(T,1);
ynew(1:K,1)=yk;

%% tracking in Fourier domain
dt=0.005;
beta=1;
for t=K:T-1

%     test=real(phi_tilde*wt);
%     figure;hold on;
%     plot(test);
%     plot(yk(1:t));
    

    sigma=GMModel.Sigma(:,:,jindex);
    mu=GMModel.mu(jindex,:);
    mu_tilde=mu(1)+1i*mu(2);
    sigmainv=diag(inv(sigma));
    sigmainv_tilde=sigmainv(1)+1i*sigmainv(2);


    wt_next=wt+dt*beta*sigmainv_tilde.*(w_tilde-wt);
    %w_traj(t-K+1)=wt_next;

    phi_next_tilde=phi_tilde(t+1,:);
    ynew(t+1)=real(phi_next_tilde*wt_next);


    wt=wt_next;
end

figure; hold on;
x=linspace(1,T,T)';
plot(x,y,'r')
plot(x,ynew,'b')
