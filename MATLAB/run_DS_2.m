%% Setting Data
close all
clear



%%
modulo_print=250;
% load('datatest12.mat', 'traj')  %st=3
% load('datatest13.mat', 'traj') %st=3
load('ImpossibleDataset_EightShape.mat', 'traj') %st=9
initial_parameters = [];
initial_parameters.M = 300;
initial_parameters.R = 1;
dirsgn = 1;
no_gauss = 20;
% data files in 3D with 'time' as first column
% Select the maximum length of the trajectory
len = 6000;
% Select how many samples to remove at beginning
st = 9;
rat = ceil(max(len/traj.size, 1));
yfar1 = interp(traj.data(st:end-st,1),rat);
yfar2 = interp(traj.data(st:end-st,2),rat);
tfar1= interp(traj.time(st:end-st,1),rat);
if size(yfar1,1) >= len
    yfar1 = yfar1(1:len,:);
    yfar2 = yfar2(1:len,:);
    tfar1 = tfar1(1:len,:);
end
tfar1 = linspace(tfar1(1),tfar1(end),size(tfar1,1));
%%
% Now we prepare the data and convert it to spherical coordinates:


type = 3;
X = [yfar1,yfar2];
smoothing = 0;
time = tfar1;
if (isfield(traj,'size'))
    T = size(yfar1,1);
else
    T = [];
end

% Get data with preferred structure
[Xdata,Xvel,Rdata,Rvel,dt,T,N,m,begin] = prepareData(type,X,time,smoothing,T);
%% Plot data
%%
plotData(Xdata,Xvel,Rdata,Rvel,T,m,'cartp');
legend('Data trajectory','Location','northeast');
%% Find $\omega$ limit set
% Initialize and perform Expectation Maximization on a Gaussian Mixture Model
% with 1 or 2 Gaussians to find the model of the limit set:
%%
j=1;
[Priors, Mu, Sigma] = EM_init_kmeans([Xdata';Xvel'], j);
[Priors, Mu, Sigma] = EM([Xdata';Xvel'], Priors, Mu, Sigma);

%% If wanted, perfom dimensionality reduction
% Decide if you want to perform dimensionality reduction:
%%
% prompt = 'Enter 1 to perform dimensionality reduction and 0 to skip it:\n';
% dimred = input(prompt);
dimred = 1;
%%
% If wanted, perform dimensionality reduction (PCA) and plot projected data:

if (j == 2)     % if there are more clusters, choose which one to perform eigenvalue decomposition on
    prompt = 'Choose which gaussian to keep (1 for 1st, 2 for 2nd):\n';
    k = input(prompt);
else
    k = 1;
end
if (dimred == 1)
    % Get rotation matrix from PCA performed on covariance matrix of gaussian:
    [Rrot,~] = eig(Sigma(1:N,1:N,k));
    Rrot = Rrot(:,N:-1:1);

    % Plot projected (rotated) data -- takes a lot of time! If you trust it, set "if false"
    if false
        figure; hold on; grid on;
        subplot(1,2,1); hold on; grid on;
        title('Original data'); xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
        subplot(1,2,2); hold on; grid on;
        title('Projected data'); xlabel('e_1'); ylabel('e_2'); zlabel('e_3');

        for i = 1:sum(T)
            Xplot = (Rrot \ (Xdata(i,:)' - Mu(1:N,k)))';
            if N == 2
                subplot(1,2,1); plot(Xdata(i,1), Xdata(i,2), 'r.'); grid on; hold on;
                subplot(1,2,2); plot(Xplot(1), Xplot(2), 'r.'); grid on; hold on;
            else
                subplot(1,2,1); view(3); plot3(Xdata(i,1), Xdata(i,2), Xdata(i,3), 'r.'); grid on; hold on;
                subplot(1,2,2); view(3); plot3(Xplot(1), Xplot(2), Xplot(3), 'r.'); grid on; hold on;
            end
        end
    end
else
    % If you are not performing dimensionality reduction
    Rrot = eye(N);
end

% Find Euler angle (or rotation angle if N = 2)
if(N == 3)
    theta0 = rotm2eul(Rrot)
elseif(N == 2)
    theta0 = acos(Rrot(1,1))
end

% Get rotated data and save original data
Xdata_ = Xdata;
%% Optimization
% Select initial values of paramaters and optimize:
%%
Xdata = Xdata_;
Xdata = (Rrot \ (Xdata' - Mu(1:N,k)))';


% initial_parameters = [];
% Set initial rho0 to variance of Gaussian model
initial_parameters.rho0 = 3*mean(diag(Sigma(1:N,1:N,1)));
% Set initial values of other parameters
initial_parameters.x0 = [0 0];
% initial_parameters.M = 20;
% Select which parameters to optimize in first objective (M, rho0, (always 0 for R,) a, x0):
initial_parameters.first = [0 1 0 1 1];
% Select which parameters to optimize in second objective (M, (always 0 for rho0,) R, a, x0):
initial_parameters.second = [0 0 1 0 0];        % e.g. [1 1 1 1 0] ecludes x0

%%%%%% OPTIMIZATION FUNCTION: %%%%%%
[params] = optimizePars(initial_parameters,Xdata,dt,begin,1);
%% Get parameters learned
%%
rho0 = params.rho0;
M = params.M;
R = params.R;
a = params.a;
% Add the mean of the Gaussian model to x0
if isfield(initial_parameters,'x0') && isempty(initial_parameters.x0)
    x0 = -Mu(1:N,k)';
else
    x0 = (Rrot * params.x0' - Mu(1:N,k))';
end
params.x0 = params.x0+x0;
params.Rrot = Rrot;
disp(params);




%% Prior computation
% rd=params.rho0;
% wd=params.R;
% for i=1:1:size(Xdata_,1)
%     %r_data(i,:) = cart2hyper(a.*(params.Rrot\(Xdata_(i,:)+params.x0)')');
%     %r_dot(i)=-sqrt(params.M.*2) .* (r_data(i,1) - rd);
%     %w(i)=params.R.* exp(-4.*params.M^2.*(r_dot(i)).^2);
%     %xd_pol=[rd,r_data(i,2)];
%     %xd(i,:)= (params.Rrot*(hyper2cart(xd_pol)./a)')' - params.x0;
%
%     [r_data(i,:),dr_data(i,:)] = DS(Xdata_(i,:),params);
%     rd_data(i,:) = r_data(i,:) + dr_data(i,:)*dt;
%     Xdata_desired(i,:) = (params.Rrot*(hyper2cart(rd_data(i,:))./params.a)')' - params.x0;
%     dXdata_desired(i,:)= Xdata_desired(i,:) - Xdata_(i,:);
% end
% means=mean(Xdata_desired);
% M=[Xdata_desired,dXdata_desired];



%% process
t=(1:size(Xdata_,1))';
dt2=0.01;
VXdata_=diff(Xdata_)/dt2;
VXdata_=[VXdata_;VXdata_(size(VXdata_,1),:)];

VXdata_= imgaussfilt(lowpass(VXdata_,1e-10));


% figure; hold on;
% subplot(2,1,1); hold on;
% plot(linspace(1,size(VXdata_,1),size(VXdata_,1)),VXdata_nf(:,1),'black');
% plot(linspace(1,size(VXdata_,1),size(VXdata_,1)),VXdata_(:,1),'red');
% subplot(2,1,2); hold on;
% plot(linspace(1,size(VXdata_,1),size(VXdata_,1)),VXdata_nf(:,2),'black');
% plot(linspace(1,size(VXdata_,1),size(VXdata_,1)),VXdata_(:,2),'red');
Data=[Xdata_,VXdata_];
pos_Y=[Xdata_(2:end,:);2*Xdata_(end,:)-Xdata_(end-1,:)];
Y=[pos_Y,VXdata_];
Y=Data;

param.mean=mean(Data);
param.sigma=sqrt(var(Data));

sigma_n=(1e-2)*eye(size(Data,1));
G = Kernel_Mult(Data,Data,param)+sigma_n;
identity=eye(size(G));
G_inv=identity/G;
%mat_=G*G_inv;

pos_init=Data(1,1:2);
vel_init=Data(1,3:4);

%% Test desired
% for i =1:size(Data,1)
%     if(mod(i,100)==0)
%         disp(i)
%     end
%     K = Kernel_Mult(Data(i,:),Data,param);
%     [r,dr] = DS(Data(i,1:2),params);
%     rd = r + dr*dt;
%     X_desired = (params.Rrot*(hyper2cart(rd)./params.a)')' - params.x0;
%     dX_desired= X_desired - Data(i,1:2);
%     m=[X_desired,dX_desired];
%     desired(i,:)=m+(K*G_inv)*(Y-m);
% end
% figure; hold on;
% scatter(Data(:,1),Data(:,2),10,"filled",'black');
% scatter(desired(:,1),desired(:,2),10,"filled",'red');
% scatter(Data(:,1)+desired(:,3),Data(:,2)+desired(:,4),10,"filled",'green');

v_mean_test= mean(sqrt(Data(:,3).^2+Data(:,4).^2));
verbose_KG = true;
params.rho0=0.35; %eight_shape_dataset
params.Rrot=[-1,-0.0089;0.0089,-0.5]; %eight_shape_dataset
params.M = 2;
params.R=-0.5;
%% Loop
%eight shape
% Kp=800e-2;
% Kv=200e-3;

loopmax=4000;
Kp=800e-2;
Kv=200e-3;
figure; hold on;
scatter(Data(:,1),Data(:,2),10,"filled",'black');

n_traj=8;
norm_traj=norm(2*[pos_init+params.x0]);
for j=1:n_traj
    points(j,:)=-params.x0+[norm_traj*cos(2*pi*j/n_traj+pi/n_traj),norm_traj*sin(2*pi*j/n_traj+pi/n_traj)];
end

for j =1:n_traj
    clear state
    state(1,:)=[points(j,:),1e-10,1e-10];
    for i = 1:loopmax%[1,500,1000,1500]%
        %             pos_init=Data(i,1:2);
        %             vel_init=Data(i,3:4);
        %             state(i,:)=[pos_init,vel_init];

        if(mod(i,100)==0)
            disp(['Trajectories: ',num2str(((j-1)*loopmax+i)/(n_traj*loopmax)*100),'%'])
        end
        %state=Data(i,:);
        K = Kernel_Mult(state(i,:),Data,param);
        %     Kxx = Kernel_Mult(state(i,:),state(i,:),param);

        %PRIOR
        [r,dr] = DS(state(i,1:2),params);
        rd = r + dr*dt;
        X_desired = (params.Rrot*(hyper2cart(rd)./params.a)')' - params.x0;
        %dX_desired= v_mean_test*(X_desired - state(i,1:2))./(norm(X_desired - state(i,1:2)));
        dX_desired= X_desired - state(i,1:2);
        m=[X_desired,dX_desired];

        %         desired(i,:)=m;                    %only prior
        %     desired(i,:)=(K*G_inv)*(Y);           %without prior

        desired(i,:)=m;%+(K*G_inv)*(Y-m);           %full
        %     sigma_d(i,:,:)=Kxx-K*G_inv*K';

        if(verbose_KG & i==4700)
            KG=K*G_inv;
            figure; hold on;
            title('Importance of each datapoint')
            scatter(state(i,1),state(i,2),40,"filled",'red');
            scatter(m(1),m(2),40,"filled",'green');
            scatter(desired(i,1),desired(i,2),40,"filled",'magenta');
            scatter3(Data(:,1),Data(:,2),KG(:),5,KG(:),'filled')
            cbar = colorbar;
            drawnow;
        end
        error(i,:)=state(i,:)-desired(i,:);

        %rule of thumb Kv<0.5Kp
        %10 and 5 e-4 works for prior only
        %dataset_12
        %     Kp=30e-4;
        %     Kv=10e-4;


        acceleration_p(i,:)=-Kp*error(i,1:2);%./param.sigma(1:2);
        acceleration_v(i,:)=-Kv*error(i,3:4);%./param.sigma(3:4);
        %delta(i,:)=acceleration_p(i,:)-acceleration_v(i,:);

        acceleration(i,:)=acceleration_p(i,:)+acceleration_v(i,:);
        new_pos=state(i,1:2)+dt2*state(i,3:4);
        new_vel=state(i,3:4)+acceleration(i,:);
        state(i+1,:)=[new_pos,new_vel];

    end
    scatter(state(:,1),state(:,2),8,"filled",'red');
    scatter(state(1,1),state(1,2),20,"filled",'blue');
end
%scatter(state(1:10,1),state(1:10,2),10,"filled",'red');
%drawnow;

clear m;
for i =1:size(Data,1)
    angle=2*pi*i/size(Data,1);
    rd = [params.rho0,angle];
    X_desired = (params.Rrot*(hyper2cart(rd)./params.a)')' - params.x0;
    m(i,:)=[X_desired];
end

figure; hold on;
scatter(Data(:,1),Data(:,2),10,"filled",'black');
scatter(m(:,1),m(:,2),8,"filled",'red');
%% Vel map
n=100;
x1=min(Data(:,1));
x2=max(Data(:,1));
y1=min(Data(:,2));
y2=max(Data(:,2));

margin=1;
xlim=[x1+margin*(x1-x2),x2+margin*(x2-x1)];
ylim=[y1+margin*(y1-y2),y2+margin*(y2-y1)];

pot_X = linspace(xlim(1),xlim(2),n)';
pot_Y = linspace(ylim(1),ylim(2),n)';
for i=1:n
    disp(['Vel map: ',num2str(i/n*100),'%'])
    for j=1:n
        state_test_pos=[pot_X(i),pot_Y(j)];
        nearest_point= dsearchn(Data(:,1:2),state_test_pos);
        vel=Data(nearest_point,1:2)-state_test_pos;
        %vel2=Data(nearest_point,3:4);
        state_test=[state_test_pos,vel];
        K = Kernel_Mult(state_test,Data,param);
        %     Kxx = Kernel_Mult(state(i,:),state(i,:),param);

        %PRIOR
        [r,dr] = DS(state_test(1:2),params);
        rd = r + dr*dt;
        X_desired = (params.Rrot*(hyper2cart(rd)./params.a)')' - params.x0;
        dX_desired= X_desired - state_test(1:2);
        m=[X_desired,dX_desired];

        %         desired=m;
        desired=m+(K*G_inv)*(Y-m);           %full
        %         error=state_test-desired;





        val=(i-1)*n+j;
        if(verbose_KG & val==1)
            KG=K*G_inv;
            KG(KG<0)=0;
            figure; hold on;
            title('Importance of each datapoint')
            scatter(state_test(1),state_test(2),40,"filled",'red');
            scatter(m(1),m(2),40,"filled",'green');
            scatter(desired(1),desired(2),40,"filled",'magenta');
            scatter3(Data(:,1),Data(:,2),KG(:),5,KG(:),'filled')
            cbar = colorbar;
            drawnow;
        end


        vel_no_acc=-(Kp*(state_test(1:2)-desired(1:2)))/(Kv)+desired(3:4);
        Ygmmbase(val,:)=[state_test(3),state_test(4)];


        %acceleration=acceleration_p+acceleration_v;
        %state_test(3:4)=state_test(3:4)+acceleration;
        %Ygmm(val,:)=[desired(3),desired(4)];
        Ygmm(val,:)=[vel_no_acc(1),vel_no_acc(2)];
        %         Ygmm(val,:)=[new_vel(1),new_vel(2)];
    end
end
figure;hold on;
title('Position desired base')
x_mesh=meshgrid(pot_X);
y_mesh=(meshgrid(pot_Y))';
scatter(Data(:,1),Data(:,2),10,"filled",'black');
streamslice(x_mesh,y_mesh,reshape(Ygmmbase(:,1),n,n),...
    reshape(Ygmmbase(:,2),n,n),2,'method','cubic');

figure;hold on;
title('Velocity desired')
x_mesh=meshgrid(pot_X);
y_mesh=(meshgrid(pot_Y))';
scatter(Data(:,1),Data(:,2),10,"filled",'black');
streamslice(x_mesh,y_mesh,reshape(Ygmm(:,1),n,n),...
    reshape(Ygmm(:,2),n,n),2,'method','cubic');

