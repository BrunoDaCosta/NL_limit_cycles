close all
clear all
load('ImpossibleDataset.mat', 'traj');
% load('datatest12.mat', 'traj');
 
data = traj.data(5:length(traj.data),:);
mean1=mean(data(:,1));
mean2=mean(data(:,2));
data3=[data(:,1).^2,data(:,2).^2,data(:,1).*data(:,2)];
%figure;
%plot3(data3(:,1),data3(:,2),data3(:,3));

kernel = Kernel('type', 'gaussian', 'gamma', 1000);
parameter = struct('numComponents', 3, 'kernelFunc', kernel);
% build a KPCA object
kpca = KernelPCA(parameter);
% train KPCA model
kpca.train(data3);

%　mapping data
mappingData = kpca.score;

% Visualization
kplot = KernelPCAVisualization();
% visulize the mapping data
kplot.score(kpca)