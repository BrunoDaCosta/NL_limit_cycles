%% Setting Data
close all
clear

plots=true;
load('datatest12.mat', 'traj');
erased_values=20;
data=traj.data(erased_values:length(traj.data),:);
X_=data(:,1);
meanX=mean(X_);
Y_=data(:,2);
meanY=mean(Y_);

X=X_-meanX;
Y=Y_-meanY;
[THETA,RHO]=cart2pol(X,Y);

%% Calculate Omega
t=(1:length(THETA))';
dt=ones(length(THETA)-1,1);
dTHETA=diff(THETA);
OMEGA=dTHETA./dt;
OMEGA=[OMEGA;OMEGA(length(OMEGA))];
while(max(OMEGA)>pi)
    OMEGA(OMEGA>pi)=OMEGA(OMEGA>pi)-2*pi;
end
while(min(OMEGA)<-pi)
    OMEGA(OMEGA<-pi)=OMEGA(OMEGA<-pi)+2*pi;
end

%% Filter (smoothing) Omega
sigma = 5;
sz = 11;    % length of gaussFilter vector
line = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-line .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter);
OMEGA= filter (gaussFilter,1, OMEGA);

%% Adding data for THETA<pi and THETA >pi
perc_b=-pi+0.52;
perc_a=pi-0.52; %(30 last degrees)
THETA_new=[THETA(THETA>perc_a)-2*pi;THETA(THETA<perc_b)+2*pi];
RHO_new=[RHO(THETA>perc_a);RHO(THETA<perc_b)];
OMEGA_new=[OMEGA(THETA>perc_a);OMEGA(THETA<perc_b)];
THETA=[THETA;THETA_new];
RHO=[RHO;RHO_new];
OMEGA=[OMEGA;OMEGA_new];

%% Precomputing GPR constant matrices

kfcn = @(XN,XM,param) param(1)*exp(-(pdist2(XN,XM).^2)/(2*param(2)^2));
param=[1,0.2];

mdl_rho = fitrgp(THETA,RHO,'KernelFunction',kfcn,'KernelParameters',param);
param_rho.g=mdl_rho.KernelInformation.KernelParameters(1);
param_rho.l=mdl_rho.KernelInformation.KernelParameters(2);
param_rho.sigma = mdl_rho.Sigma;
disp("optimization of hyperparameters for E(Rho|Theta) DONE!");

mdl_omega = fitrgp([THETA,RHO],OMEGA,'KernelFunction',kfcn,'KernelParameters',param);
param_omega.g =mdl_omega.KernelInformation.KernelParameters(1);
param_omega.l =mdl_omega.KernelInformation.KernelParameters(2);
param_omega.sigma = mdl_omega.Sigma;
disp("optimization of hyperparameters for E(Omega|Rho,Theta) DONE!");



Bg = Gauss_Regression_Precompute(THETA,param_rho);
Bh = Gauss_Regression_Precompute([RHO,THETA],param_omega);

%% Applying hand-made GPR
rho_pred = Gauss_Regression(Bg,RHO,THETA,THETA,param_rho);
ome_pred = Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO,THETA],param_omega);

%% Calculating the limit cycle
steps=500;
RHO_vect=zeros(steps,1);
THETA_vect=zeros(steps,1);
OMEGA_vect=zeros(steps,1);
X_vect=zeros(steps,1);
Y_vect=zeros(steps,1);
RHO_vect(1,1)=RHO(1,1);
THETA_vect(1,1)=THETA(1,1);

X_vect(1,1)=1.5*X(1);
Y_vect(1,1)=1.5*Y(1);
dX=zeros((steps)-1,1);
dY=zeros((steps)-1,1);
for i = 1:(steps)-1
    [dRHO,dOMEGA] = next_step(RHO,THETA,OMEGA,RHO_vect(i,1),THETA_vect(i,1),OMEGA_vect(i,1),Bg,Bh,param_rho,param_omega);
    
    %RHO_vect(i+1,1)=RHO_vect(i,1) + dRHO;
    OMEGA_vect(i+1,1)=OMEGA_vect(i,1) + dOMEGA;
    %THETA_vect(i+1,1)=THETA(i,1) + OMEGA_vect(i+1,1);
    
    dX(i)=dRHO*cos(THETA_vect(i,1))-RHO_vect(i,1)*OMEGA_vect(i+1,1)*sin(THETA_vect(i,1));
    dY(i)=dRHO*sin(THETA_vect(i,1))+RHO_vect(i,1)*OMEGA_vect(i+1,1)*cos(THETA_vect(i,1));
    
    X_vect(i+1,1)=X_vect(i,1)+dX(i);
    Y_vect(i+1,1)=Y_vect(i,1)+dY(i);
    [THETA_vect(i+1,1),RHO_vect(i+1,1)]=cart2pol(X_vect(i+1,1),Y_vect(i+1,1));
end
disp("Limit cycle calculated");

%% Testing/Plotting
if(plots)
    figure;
    scatter(X,Y,5,'filled','r');
    %
    %     figure;
    %     scatter(THETA,RHO,5,'filled');
    
    figure; hold on;
    scatter(THETA,rho_pred,5,'filled','g');
    scatter(THETA,RHO,5,'filled','r');
    legend({'rho predicted','rho'})
    xlabel('THETA');
    ylabel('RHO');
    
    figure;
    scatter3(RHO,THETA,ome_pred,'filled','g');
    hold on;
    scatter3(RHO,THETA,OMEGA,'filled','r');
    legend({'w predicted','w'})
    xlabel('RHO');
    ylabel('THETA');
    zlabel('OMEGA');
    hold off;
    
    figure; hold on;
    scatter(X,Y,5,'filled','r');
    scatter(X_vect,Y_vect,20,'filled','black');
    legend({'Dataset','Path limit cycle'})
    xlabel('X');
    ylabel('Y');
end

%% Create potential map
n=200;
X1 = min(X)+0.05*(min(X)-max(X));
X2 = max(X)+0.05*(max(X)-min(X));
Y1 = min(Y)+0.05*(min(Y)-max(Y));
Y2 = max(Y)+0.05*(max(Y)-min(Y));
pot_X = linspace(X1,X2,n)';
pot_Y = linspace(Y1,Y2,n)';
pot_Z = zeros(n,n);
N=1;
T=1;
K_alpha=5;
for i=1:n
    disp(i/n);
    for j=1:n
        [THETA_,RHO_]=cart2pol(pot_X(i),pot_Y(j));
        OMEGA_desired=Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO_,THETA_],param_omega);
        RHO_desired=Gauss_Regression(Bg,RHO,THETA,THETA_,param_rho);
        
        alpha=min(K_alpha*(abs(RHO_-RHO_desired)),1);
        OMEGA_= (1-alpha)*OMEGA_desired;
        
        
        
        [dRHO_,dOMEGA_] = next_step(RHO,THETA,OMEGA,RHO_,THETA_,OMEGA_,Bg,Bh,param_rho,param_omega);
        OMEGA_=OMEGA_ + dOMEGA_;
        %dX_=dRHO_*cos(THETA_)-RHO_*OMEGA_*sin(THETA_);
        %dY_=dRHO_*sin(THETA_)+RHO_*OMEGA_*cos(THETA_);
        pot_Z(j,i)=N*(dRHO_^2+RHO_^2*OMEGA_^2)+T*abs(RHO_-RHO_desired)^2;
        
    end
end


axis=[linspace(0,0.03,10),linspace(0.0333,max(pot_Z,[],'all'),10)];
figure; hold on;
pot_Z_new=pot_Z;
pot_Z_new(pot_Z_new>0.001)=0.001;
contourf(pot_X,pot_Y,pot_Z_new,30);
scatter(X_vect,Y_vect,15,'filled','r')
cbh = colorbar ; %Create Colorbar



