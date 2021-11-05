clear; close all;

%% 2D test
x=linspace(-2*pi,2*pi,1000)';
x_test=linspace(-2*pi,2*pi,74)';
y=sin(x)+0.05*rand(size(x));

%Bx = Gauss_Regression_Precompute(x);
Bx = Gauss_Regression_Precompute(x);
y_pred = Gauss_Regression(Bx,y,x,x_test);

figure;hold on;
plot(x,y,'r');
scatter(x_test,y_pred,10,'g','filled');
xlabel('x');
ylabel('y');

%% 3D test
big_data=100;
test_data=22;
x=repmat(linspace(-2*pi,2*pi,big_data),big_data,1)';
x=reshape(x,big_data^2,1);
y=repmat(linspace(0,10,big_data),big_data,1);
y=reshape(y,big_data^2,1);
z=sin(x)+sin(y).^2;

x_=repmat(linspace(-2*pi,2*pi,test_data),test_data,1)';
y_=repmat(linspace(0,10,test_data),test_data,1);

Bxy = Gauss_Regression_Precompute([x,y]);
z_pred=zeros(22);
z_pred_dx=zeros(22);
z_pred_dy=zeros(22);
for i =1:22
    z_pred(:,i) = Gauss_Regression(Bxy,z,[x,y],[x_(:,i),y_(:,i)]);
    z_pred_dx(:,i) = Gauss_Regression(Bxy,z,[x,y],[x_(:,i),y_(:,i)],1);
    z_pred_dy(:,i) = Gauss_Regression(Bxy,z,[x,y],[x_(:,i),y_(:,i)],2);
end

figure;
scatter3(x,y,z,5,'r','filled');
hold on;
scatter3(reshape(x_,test_data^2,1),reshape(y_,test_data^2,1),reshape(z_pred,test_data^2,1),10,'g','filled');
title('Predicted Z');
xlabel('x');
ylabel('y');
zlabel('z');

figure;
scatter3(x,y,z,5,'r','filled');
hold on;
scatter3(reshape(x_,test_data^2,1),reshape(y_,test_data^2,1),reshape(z_pred_dx,test_data^2,1),10,'g','filled');
title('Predicted dZ/dx');
xlabel('x');
ylabel('y');
zlabel('z');

figure;
scatter3(x,y,z,5,'r','filled');
hold on;
scatter3(reshape(x_,test_data^2,1),reshape(y_,test_data^2,1),reshape(z_pred_dy,test_data^2,1),10,'g','filled');
title('Predicted dZ/dy');
xlabel('x');
ylabel('y');
zlabel('z');
