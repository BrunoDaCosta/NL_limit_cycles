%%
clear all;
close all;

load('datatest12.mat', 'traj')
initial_parameters = [];
initial_parameters.M = 300;
initial_parameters.R = 1; 
dirsgn = 1;
no_gauss = 20;
% data files in 3D with 'time' as first column
    % Select the maximum length of the trajectory
    len = 6000;
    % Select how many samples to remove at beginning
    st = 5;
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

%% Find $\omega$ limit set
% Initialize and perform Expectation Maximization on a Gaussian Mixture Model 
% with 1 or 2 Gaussians to find the model of the limit set:
j=1;
[Priors, Mu, Sigma] = EM_init_kmeans([Xdata';Xvel'], j);
[Priors, Mu, Sigma] = EM([Xdata';Xvel'], Priors, Mu, Sigma);

%% If wanted, perfom dimensionality reduction
% Decide if you want to perform dimensionality reduction:
dimred = 1;

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
params.x0 = x0;
params.Rrot = Rrot;
disp(params);

%% Nonlinear Limit Cycle
ZData = Xdata;

thData = zeros(size(ZData,1),1);
rsData = zeros(size(ZData,1),1);
for cnt=1:size(ZData,1)
    rr = cart2hyper(ZData(cnt,:).*a);
    [theta,~]= cart2pol(ZData(cnt,1)*a(1),ZData(cnt,2)*a(2));
    thData(cnt) = theta;
    r0 = rho0;
    rsData(cnt) = rr(:,1) / r0;
end

[t_sort,t_order] = sort(thData);

rdata = [t_sort,rsData(t_order,1)];
% ----
window_size = 60;
bins = linspace( min(t_sort),max(t_sort), window_size);
mu_bin = zeros(length(bins)-1,2);
sig_bin = zeros(length(bins)-1,2);
count_bin = zeros(length(bins)-1,1);
for ic = 1:length(bins)-1
    logicalindex = (t_sort >bins(ic)) & (t_sort < bins(ic+1));
    indx = find(logicalindex);
    bin_data = [t_sort(indx,1),rdata(indx,2)];
    bin_slice(ic).data = bin_data;
    count_bin(ic) = length(bin_data);
    mu_bin(ic,:) = mean(bin_data);
    sig_bin(ic,:) = std(bin_data);
end
for iter=1:10
    meancount = ceil(.9*mean(count_bin));
    for jc = 1:length(bin_slice)
        if (count_bin(jc) < meancount)
            sigma_ = .1*[min(sig_bin(jc,1),.01) 0;0 min(sig_bin(jc,2),.01)];
            newpoints = mvnrnd(mu_bin(jc,:),sigma_*sigma_,meancount - count_bin(jc));
            bin_slice(jc).data = [bin_slice(jc).data;newpoints];
            count_bin(jc) = count_bin(jc) + length(newpoints);
            mu_bin(jc,:) = mean(bin_slice(jc).data);
            sig_bin(jc,:) = std(bin_slice(jc).data);
        end
    end
end

rsdata = [];
for kc = 1:length(bin_slice)
    rsdata =[rsdata;bin_slice(kc).data];
end

[tt_sort,tt_order] = sort(rsdata(:,1));
rrdata = [tt_sort,rsdata(tt_order,2)];
% ---

addfirst = (rrdata(end-500:end,:)'+ [rrdata(1,1);0]-[rrdata(end,1);0])';
addlast = (rrdata(1:500,:)'- [rrdata(1,1);0]+[rrdata(end,1);0])';
Rdata = [addfirst;rrdata ; addlast];

options_gmm = statset('Display','final','MaxIter',1500,'TolFun',1e-8);
rho_GMModel = fitgmdist(Rdata,12,'Options',options_gmm);

%% Plot learned dynamics with original data and new trajectories

% Get unrotated data
Xdata = Xdata_;

% Plot original data
figure; hold on; grid on;
% plot(Xdata(:,1),Xdata(:,2),'r.'); hold on;


% Test dynamics for T time steps
X0 = 1.5*Xdata(1,:);
time = ceil((3/4)*T(1));
speed=zeros(time,1);
speedgmm=zeros(time,1);
Xplot=zeros(time,2);
Xplotgmm=zeros(time,2);
for j = 1:size(X0,1)
    X = X0(j,:);
    Xgmm = X0(j,:);
    for i = 1:time
        %%%%%% FUNCTION TO CALCULATE POLAR/SPHERICAL VELOCITIES: %%%%%%
        [r,dr] = DS(X,params);
        [rgmm,drgmm] = DSGMM(Xgmm,params,rho_GMModel);

        %%%%%% INTEGRATE THE DS TO GET NEXT POLAR/SPHERICAL POSITION: %%%%%%
        next_r = r + dr*dt;
        next_rgmm = rgmm + drgmm*dt;
        % Get next cartesian position
        X = (Rrot*(hyper2cart(next_r)./a)')' - x0;
        Xgmm = (Rrot*(hyper2cart(next_rgmm)./a)')' - x0;

        Xplot(i,:)=X;
        Xplotgmm(i,:)=Xgmm;
    end
end

first_values=50;

for i = 2:time
    speed(i-1)=sqrt((Xplot(i,1)-Xplot(i-1,1))^2+(Xplot(i,2)-Xplot(i-1,2))^2);
    speedgmm(i-1)=sqrt((Xplotgmm(i,1)-Xplotgmm(i-1,1))^2+(Xplotgmm(i,2)-Xplotgmm(i-1,2))^2);
end
plot(Xplot(:,1),Xplot(:,2),'r','MarkerSize',10); hold on; grid on;
plot(Xplotgmm(:,1),Xplotgmm(:,2),'b','MarkerSize',10); hold on; grid on;
legend('Position linear', 'Non-Linear')


figure; hold on; grid on;
plot(speed([first_values:size(speed)]),'r')
plot(speedgmm([first_values:size(speedgmm)]),'b')
legend('velocity linear', 'velocity GMM')
xlabel('Timesteps') 
ylabel('Meters / Timestep') 

%%

if(true) %plot the countourf for velocitites
    xmin_=min(Xplot(:,1));
    xmax_=max(Xplot(:,1));
    ymin_=min(Xplot(:,2));
    ymax_=max(Xplot(:,2));
    xmin=xmin_+0.1*(xmin_-xmax_);
    xmax=xmax_+0.1*(xmax_-xmin_);
    ymin=ymin_+0.1*(ymin_-ymax_);
    ymax=ymax_+0.1*(ymax_-ymin_);
    nmb=100;
    [Xmesh,Ymesh] = meshgrid(xmin:(xmax-xmin)/nmb:xmax,ymin:(ymax-ymin)/nmb:ymax);
    Zmeshlin=zeros(nmb+1);
    Zmeshnonlin=zeros(nmb+1);
    for i=1:nmb+1
        disp(i)
        for j=1:nmb+1
            [~,val]=DS([Xmesh(1,i),Ymesh(j,1)],params);
            val=(Rrot*(hyper2cart(val)./a)')' - x0;
            Zmeshlin(i,j)=sqrt(val(1)^2+val(2)^2);

            [~,val]=DSGMM([Xmesh(1,i),Ymesh(j,1)],params,rho_GMModel);
            val=(Rrot*(hyper2cart(val)./a)')' - x0;
            Zmeshnonlin(i,j)=sqrt(val(1)^2+val(2)^2);
        end
    end
    figure; hold on; grid on;
    contourf(Xmesh,Ymesh,Zmeshlin,20)
    ylim([ymin,ymax]);
    xlim([xmin,xmax]);
    
    figure; hold on; grid on;
    contourf(Xmesh,Ymesh,Zmeshnonlin,20)
    ylim([ymin,ymax]);
    xlim([xmin,xmax]);
end



