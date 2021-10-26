clear; close all;

%% 2D test
x=linspace(-2*pi,2*pi,1000)';
x_test=linspace(-2*pi,2*pi,74)'
y=sin(x)+0.05*rand(size(x));

Bx = Gauss_Regression_Precompute(x);
y_pred = Gauss_Regression(Bx,y,x,x_test);

figure;hold on;
plot(x,y,'r');
scatter(x_test,y_pred,10,'g','filled');

%% 3D test
x=linspace(-2*pi,2*pi,1000)';
y=linspace(-2*pi,2*pi,1000)';
z=sin(x)+sin(y);


Bxy = Gauss_Regression_Precompute([x,y]);
z_pred = Gauss_Regression(Bxy,z,[x,y],[x,y]);

figure;
plot3(x,y,z,'r');
hold on;
scatter3(x,y,z_pred,10,'g','filled');

