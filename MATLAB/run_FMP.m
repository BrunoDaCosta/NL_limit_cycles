clear; close all;

%% Creation of dataset
plots=false;


K=120;
load('datatest12.mat', 'traj');
period=346;
first_values_erased=3;
T=period;
Length_data=size(traj.data,1);
x_=traj.data(first_values_erased:end,1);
y_=traj.data(first_values_erased:end,2);
extension=5;
Textended=extension*T;

% Same initial values or another one
xk=(2*max(abs(x_-mean(x_))))*cos(0.1*linspace(1,K,K)')+mean(x_);
yk=(2*max(abs(y_-mean(y_))))*-sin(0.1*linspace(1,K,K)')+mean(y_);

%% partial trajectory mapping
%taking all dataset
% for i=0:2*T+1
%     begin=i*floor(Length_data/((4)*T+1));
%     x(:,i+1)=[x_(begin+1:begin+T)];
%     y(:,i+1)=[y_(begin+1:begin+T)];
% end

%taking a part repeatedy
for i=1:4*T+1
    begin=max(mod(i,T+1),1);
    shift=mod(i,5);
    x(:,i)=[x_(1+shift:T+shift)];%-0.005+0.01*rand(size(x_(1:T)));
    y(:,i)=[y_(1+shift:T+shift)];%-0.005+0.01*rand(size(y_(1:T)));
%     x(:,i)=[x_(begin+1:T);x_(1:begin)];
%     y(:,i)=[y_(begin+1:T);y_(1:begin)];
end
% plot to see 10 different N
% plot(x(:,1:100:1000))

%data_pos=[x;y];

[E_mu_x,w_data_x,muj_x,phi_x,wt_x,jindex_x,GMModel_x]=partial_trajectory_mapping(x,xk,plots);
xnew=zeros(Textended,1);
xnew(1:K,1)=xk;

[E_mu_y,w_data_y,muj_y,phi_y,wt_y,jindex_y,GMModel_y]=partial_trajectory_mapping(y,yk,plots);
ynew=zeros(Textended,1);
ynew(1:K,1)=yk;





%% tracking in Fourier domain
dt=0.005;
beta=0.05;
alpha=0.2;
saturation_x=mean(max(abs(diff(x))));
saturation_y=mean(max(abs(diff(y))));
t=K;
sigmainv_x=eye(size(GMModel_x.Sigma(:,:,jindex_x)))/GMModel_x.Sigma(:,:,jindex_x);
sigmainv_y=eye(size(GMModel_y.Sigma(:,:,jindex_y)))/GMModel_y.Sigma(:,:,jindex_y);
sigmainv_x=diag(diag(sigmainv_x));
sigmainv_y=diag(diag(sigmainv_y));

for ttotal=K:Textended-1
    if(mod(ttotal,10)==0)
        disp(ttotal);
    end

    if(t+1>=T)
        t=0;

%         [E_mu_x,w_data_x,muj_x,phi_x,wt_x,jindex_x,GMModel_x]=partial_trajectory_mapping(x,xnew(ttotal-K+1:ttotal),plots);
%         [E_mu_y,w_data_y,muj_y,phi_y,wt_y,jindex_y,GMModel_y]=partial_trajectory_mapping(y,ynew(ttotal-K+1:ttotal),plots);
%         sigmainv_x=eye(size(GMModel_x.Sigma(:,:,jindex_x)))/GMModel_x.Sigma(:,:,jindex_x);
%         sigmainv_y=eye(size(GMModel_y.Sigma(:,:,jindex_y)))/GMModel_y.Sigma(:,:,jindex_y);
%         sigmainv_x=diag(diag(sigmainv_x));
%         sigmainv_y=diag(diag(sigmainv_y));
        muj_x=E_mu_x;
        muj_y=E_mu_y;
            
    end





    %next position x
    wt_next_x=wt_x+dt*beta*sigmainv_x*(muj_x-wt_x);
    phi_next_x=phi_x(t+1,:);
    xnew(ttotal+1)=alpha*real(phi_next_x*wt_next_x)+(1-alpha)*xnew(ttotal);

    %saturation x
    if(xnew(ttotal+1)-xnew(ttotal)>saturation_x)
        xnew(ttotal+1)=xnew(ttotal)+saturation_x;
    elseif (xnew(ttotal+1)-xnew(ttotal)<-saturation_x)
        xnew(ttotal+1)=xnew(ttotal)-saturation_x;
    end

    

    %next position y
    wt_next_y=wt_y+dt*beta*sigmainv_y*(muj_y-wt_y);
    phi_next_y=phi_y(t+1,:);
    ynew(ttotal+1)=alpha*real(phi_next_y*wt_next_y)+(1-alpha)*ynew(ttotal);

    %saturation y
    if(ynew(ttotal+1)-ynew(ttotal)>saturation_y)
        ynew(ttotal+1)=ynew(ttotal)+saturation_y;
    elseif (ynew(ttotal+1)-ynew(ttotal)<-saturation_y)
        ynew(ttotal+1)=ynew(ttotal)-saturation_y;
    end

    delta_x(ttotal+1-K)=norm(muj_x-wt_x);
    delta_y(ttotal+1-K)=norm(muj_y-wt_y);

    wt_x=wt_next_x;
    wt_y=wt_next_y;

    t=t+1;
end



% truc=real(phi_x*muj_x);
% figure; hold on;
% plot(truc);
% plot(x(:,1))
% legend({'Phi*Mu_j','First demonstration'})
% title('phi*muj')


figure;
%x plot
subplot(2,2,1);hold on;
plot(x(:,1),'r')
plot(xnew,'b')
legend({'Dataset','Computed X'}, 'Location','southeast')
title('X axis')
hold off;

%y plot
subplot(2,2,2);hold on;
plot(y(:,1),'r')
plot(ynew,'b')
legend({'Dataset','Computed Y'}, 'Location','southeast')
title('Y axis')
hold off;

%full plot
subplot(2,2,[3,4]);hold on;
scatter(x(:,1),y(:,1),10,'filled','r')
plot(xnew,ynew,'b')
legend({'Dataset','Computed position'}, 'Location','northwest')
title('Robot movement')
hold off;

%delta along time
figure; hold on;
title('||omega_t - omega_d||^2')
plot(delta_x);
plot(delta_y);
legend({'X-axis','Y-axis'})


%%Presentation plots
figure; hold on; 
plot(xnew(1:K),'b'); 
plot(x(:,1),'r');
legend({'initial movement','dataset'})

figure; hold on; 
plot(xnew(1:K),'b'); 
plot(x(:,1),'r');
shift=K:size(x,1);
plot(shift,xnew(K:size(x,1)),'black'); 
legend({'initial movement','dataset','computed movement'})




