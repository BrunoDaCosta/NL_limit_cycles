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

%% Adding data for THETA<pi and THETA >pi
full_cycle=false;
negative=false;
theta_init=THETA(1);
i=0;
while(not(full_cycle))
    i=i+1;
    if(negative)
        if(THETA(i)>theta_init)
            full_cycle=true;
            i=i-1;
        end
    elseif(THETA(i)<0)
        negative=true;
    end
end

THETA=[THETA(1:i)-2*pi;THETA(1:i)+2*pi;THETA];
RHO=[RHO(1:i);RHO(1:i);RHO];

%% Preprocess Omega
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

sigma = 5;
sz = 11;    % length of gaussFilter vector
line = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-line .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter);
OMEGA= filter (gaussFilter,1, OMEGA);

%% Precomputing GPR constant matrices
Bg = Gauss_Regression_Precompute(THETA);
Bh = Gauss_Regression_Precompute([RHO,THETA]);
Bh_rho = Bh(1:size(Bh,1)/2,:);
Bh_theta = Bh(size(Bh,1)/2+1:size(Bh,1),:);

%% Applying hand-made GPR
rho_pred = Gauss_Regression(Bg,RHO,THETA,THETA);
ome_pred = Gauss_Regression(Bh,OMEGA,[RHO,THETA],[RHO,THETA]);

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
    [dRHO,dOMEGA] = next_step(RHO,THETA,OMEGA,RHO_vect(i,1),THETA_vect(i,1),OMEGA_vect(i,1),Bg,Bh,Bh_rho,Bh_theta);
    
    %RHO_vect(i+1,1)=RHO_vect(i,1) + dRHO;
    OMEGA_vect(i+1,1)=OMEGA_vect(i,1) + dOMEGA;
    %THETA_vect(i+1,1)=THETA(i,1) + OMEGA_vect(i+1,1);

    dX(i)=dRHO*cos(THETA_vect(i,1))-RHO_vect(i,1)*OMEGA_vect(i+1,1)*sin(THETA_vect(i,1));
    dY(i)=dRHO*sin(THETA_vect(i,1))+RHO_vect(i,1)*OMEGA_vect(i+1,1)*cos(THETA_vect(i,1));

    X_vect(i+1,1)=X_vect(i,1)+dX(i);
    Y_vect(i+1,1)=Y_vect(i,1)+dY(i);
    [THETA_vect(i+1,1),RHO_vect(i+1,1)]=cart2pol(X_vect(i+1,1),Y_vect(i+1,1));
end

%% Testing/Plotting
if(plots)
%     figure; 
%     scatter(X,Y,5,'filled');
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
    scatter(X_vect,Y_vect,20,'filled','g');
    legend({'Dataset','Path limit cycle'})
    xlabel('X');
    ylabel('Y');
end

%% Create potential map
n=10;
X1 = min(X)+0.05*(min(X)-max(X));
X2 = max(X)+0.05*(max(X)-min(X));
Y1 = min(Y)+0.05*(min(Y)-max(Y));
Y2 = max(Y)+0.05*(max(Y)-min(Y));
pot_X = linspace(X1,X2,n);
pot_Y = linspace(Y1,Y2,n);
pot_Z = zeros(n,n);
for i=1:n
    i=i
    for j=1:n
        OMEGA_=0.02;
        [THETA_,RHO_]=cart2pol(pot_X(i),pot_Y(j));
        [dRHO_,dOMEGA_] = next_step(RHO,THETA,OMEGA,RHO_,THETA_,OMEGA_,Bg,Bh,Bh_rho,Bh_theta);
        OMEGA_=OMEGA_ + dOMEGA_;
        dX_=dRHO_*cos(THETA_)-RHO_*OMEGA_*sin(THETA_);
        dY_=dRHO_*sin(THETA_)+RHO_*OMEGA_*cos(THETA_);

        pot_Z(i,j)=sqrt(dX_^2+dY_^2);
    end
end
figure; hold on;
contourf(pot_X,pot_Y,pot_Z);
colorbar;




