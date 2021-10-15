
clear all;
close all;

FileData = load('tab1_1.mat');
vectorzero=zeros(size(FileData.traj.time));
data=[FileData.traj.time,FileData.traj.data,vectorzero];
csvwrite('tab1_1.csv', data);