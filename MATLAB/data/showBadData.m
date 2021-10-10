%%
clear all;
close all;

%%
load('ImpossibleDataset.mat', 'traj')
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

%% Plot dataset
plotData(Xdata,Xvel,Rdata,Rvel,T,m,'cartp');
legend('Data trajectory','Location','northeast');

%% Find $\omega$ limit set
% Initialize and perform Expectation Maximization on a Gaussian Mixture Model 
% with 1 or 2 Gaussians to find the model of the limit set:
j=1;
[Priors, Mu, Sigma] = EM_init_kmeans([Xdata';Xvel'], j);
[Priors, Mu, Sigma] = EM([Xdata';Xvel'], Priors, Mu, Sigma);

%% If wanted, perfom dimensionality reduction
% If wanted, perform dimensionality reduction (PCA) and plot projected data:
k=1;
% Get rotation matrix from PCA performed on covariance matrix of gaussian:
[Rrot,~] = eig(Sigma(1:N,1:N,k));
Rrot = Rrot(:,N:-1:1);
theta0 = acos(Rrot(1,1))

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
%             randindex = ceil(length(bin_slice(jc).data)*rand(meancount - count_bin(jc),1) );
%             newpoints =  bin_slice(jc).data(randindex,:);
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
gmPDF = @(x,y)reshape(pdf(rho_GMModel,[x(:) y(:)]),size(x));

figure;
scatter(Rdata(:,1), Rdata(:,2),5,'MarkerEdgeColor','none',...
    'MarkerFaceColor','r','LineWidth',.5);
hold on;
fcontour(gmPDF,[-4 4 .3 1.5],'LevelList',.1:.2:1.5,'LineWidth', 2);
grid on;
ax = gca;
ax.FontSize = 14; 
ylabel('$\bf{s}(\theta)$','Interpreter','latex','FontSize',18,'FontWeight','bold');
xlabel('$\theta$','Interpreter','latex','FontSize',18,'FontWeight','bold');
ylim([.3 1.7]);
xlim([-4. 4]);
hold off;

%% Plot learned dynamics with original data and new trajectories

% Get unrotated data
Xdata = Xdata_;

% Plot original data
figure; hold on; grid on;

plot(Xdata(:,1),Xdata(:,2),'r.'); hold on;

% Plot streamlines / arrows to show dynamics
yl = ylim;
xl = xlim;
resol = 500;
[Xs,Ys] = meshgrid(linspace(xl(1),xl(2),resol),linspace(yl(1),yl(2),resol));
X_plot = [Xs(:), Ys(:)];

% Test dynamics for T time steps
X0 = 1.5*Xdata(1,:);
colorrgb = 0;
nmbpoints=round(T(1)/1.5);
X=X0;
for i = 1:nmbpoints
    X_prev = X;
    %%%%%% FUNCTION TO CALCULATE POLAR/SPHERICAL VELOCITIES: %%%%%%
    [r,dr] = DS(X,params);
    %%%%%% INTEGRATE THE DS TO GET NEXT POLAR/SPHERICAL POSITION: %%%%%%
    next_r = r + dr*dt;
    % Get next cartesian position
    X = (Rrot*(hyper2cart(next_r)./a)')' - x0;
    if N == 2
        %plot(X(1),X(2),'k.','LineWidth',2,'MarkerSize',10,'Color',[0.0 colorrgb 0.9999-colorrgb]); hold on; grid on;
        plot(X(1),X(2),'k.','LineWidth',2,'MarkerSize',10,'Color',[0.0 0.0 0.0]); hold on; grid on;
    else
        plot3(X(1),X(2),X(3),'k.'); hold on; grid on;
    end

    %color (Blue to green)
    %colorrgb = colorrgb + 1/(nmbpoints+1);
    if(i>=40)
        if(X(1)<xl(1))
            xl(1)=X(1)
        end
        if(X(1)>xl(2))
            xl(2)=X(1)
        end
        if(X(2)<yl(1))
            yl(1)=X(2)
        end
        if(X(2)>yl(2))
            yl(2)=X(2)
        end
    end
end
ylim(yl);
xlim(xl);
%========================
% hold on;
fig2 = figure;
plot(Xdata(:,1),Xdata(:,2),'r.'); hold on;
Ygmm = zeros(size(X_plot));

streamslice(Xs,Ys,reshape(Ygmm(:,1),resol,resol),...
reshape(Ygmm(:,2),resol,resol),'method','cubic');
ax = gca;
ax.FontSize = 14; 
ylim(yl);
xlim(xl);
ylabel('$\xi_2$','Interpreter','latex','FontSize',18,'FontWeight','bold');
xlabel('$\xi_1$','Interpreter','latex','FontSize',18,'FontWeight','bold');
hold on;
% grid on;
%% Test dynamics (non-linear) for T time steps

%This sets the point as the first one of data (x1.5)
X0 = 1.5*Xdata(1,:);
nmbpoints=round(1.2*T(1));
colorrgb = 0;
X = X0;
%modulo to plot less points (process faster)
modulo=1;
for i = 1:(nmbpoints) 
    X_prev = X;
    %%%%%% FUNCTION TO CALCULATE POLAR VELOCITIES: %%%%%%
    [r,dr] = DSGMM(X,params,rho_GMModel);
    %%%%%% INTEGRATE THE DS TO GET NEXT POLAR POSITION: %%%%%%
    next_r = r + dr*dt;
    % Get next cartesian position
    X = (Rrot*(hyper2cart(next_r)./a)')' - x0;
    if (mod(i,modulo)==0)
        %plot(X(1),X(2),'k.','LineWidth',2,'MarkerSize',10,'Color',[(sin(-2*pi/3+colorrgb)+1)/2 (sin(colorrgb)+1)/2 (sin(2*pi/3+colorrgb)+1)/2]); hold on; grid on;
        plot(X(1),X(2),'k.','LineWidth',2,'MarkerSize',10,'Color',[0.0 0.0 0.0]); hold on; grid on;
    end
    %color (Rainbow)
    %colorrgb = colorrgb + (2*pi)/(nmbpoints+1);
end

