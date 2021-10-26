close all
clear
%load('datatest12.mat', 'traj');
load('ImpossibleDataset.mat', 'traj');
 

data2 = traj.data(5:length(traj.data),:);
%delta=mean(mean(data2));
delta=1;
x=data2(:,1);
y=data2(:,2);
z = delta*(ones(length(data2), 1));
z(32,1)=delta/2;
data=[x,y,z];


[TH,PHI,R] = cart2sph(data(:,1),data(:,2),data(:,3));

figure;
plot3(data(:,1),data(:,2),data(:,3))
xlabel('x','Interpreter','latex','FontSize',18,'FontWeight','bold');
ylabel('y','Interpreter','latex','FontSize',18,'FontWeight','bold');
zlabel('z','Interpreter','latex','FontSize',18,'FontWeight','bold');

if(true)%plot the 3 combinations of spherical coordinates
    figure;
    plot(TH,R)
    ylabel('$\rho$($\theta$)','Interpreter','latex','FontSize',18,'FontWeight','bold');
    xlabel('$\theta$','Interpreter','latex','FontSize',18,'FontWeight','bold');

    figure;
    plot(TH,PHI)
    ylabel('$\phi$($\theta$)','Interpreter','latex','FontSize',18,'FontWeight','bold');
    xlabel('$\theta$','Interpreter','latex','FontSize',18,'FontWeight','bold');

    figure;
    plot(PHI,R)
    ylabel('$\rho$($\phi$)','Interpreter','latex','FontSize',18,'FontWeight','bold');
    xlabel('$\phi$','Interpreter','latex','FontSize',18,'FontWeight','bold');
end

figure;
plot3(TH,PHI,R)
xlabel('$\theta$','Interpreter','latex','FontSize',18,'FontWeight','bold');
ylabel('$\phi$','Interpreter','latex','FontSize',18,'FontWeight','bold');
zlabel('$\rho$','Interpreter','latex','FontSize',18,'FontWeight','bold');




