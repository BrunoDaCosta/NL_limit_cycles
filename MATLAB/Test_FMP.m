clear; close all;

%% Creation of dataset
plots=false;

T=1000;
K=123;
x=linspace(0,20.5*pi,T)';
%y=dataX;
y=sin(x);
extension=3;
Textended=extension*T;

% Same initial values or another one
%yk=sin(x(1:K));
yk=0.2*cos(2*x(1:K));

%% partial trajectory mapping 
[w_tilde,phi_tilde,wt,jindex,GMModel]=partial_trajectory_mapping(y,yk,plots);
ynew=zeros(T,1);
ynew(1:K,1)=yk;


%% tracking in Fourier domain
dt=0.005;
beta=2000;
saturation=max(abs(diff(y)));
for t=K:extension*T-1

%     test=real(phi_tilde*wt);
%     figure;hold on;
%     plot(test);
%     plot(yk(1:t));
    

    sigma=GMModel.Sigma(:,:,jindex);
    mu=GMModel.mu(jindex,:);
    mu_tilde=mu(1)+1i*mu(2);
    sigmainv=diag(inv(sigma));
    sigmainv=inv([cov(w_tilde),0;0,cov(w_tilde)]);

    sigmainv_tilde=sigmainv(1)+1i*sigmainv(2);

    wt_next=wt+dt*beta*sigmainv_tilde.*(w_tilde-wt);
    %w_traj(t-K+1)=wt_next;

    phi_next_tilde=phi_tilde(mod(t+1,T)+1,:);
    ynew(t+1)=real(phi_next_tilde*wt_next);

    %saturation
    if(ynew(t+1)-ynew(t)>saturation)
        ynew(t+1)=ynew(t)+saturation;
    elseif (ynew(t+1)-ynew(t)<-saturation)
        ynew(t+1)=ynew(t)-saturation;
    end


    wt=wt_next;
end


xextended=linspace(0,extension*max(x),Textended)';
yextended=zeros(size(xextended));
yextended(1:size(y,1),:)=y;

figure; hold on;
plot(xextended,yextended,'r')
plot(xextended,ynew,'b')
legend({'Dataset','Robot arm'})

%scatter(real(w_tilde),imag(w_tilde),10,'filled','r')
