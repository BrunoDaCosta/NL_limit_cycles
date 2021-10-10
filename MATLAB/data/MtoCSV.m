
clear all;
close all;

FileData = load('ImpossibleDataset.mat');
vectorzero=zeros(size(FileData.traj.time));
data=[FileData.traj.time,FileData.traj.data,vectorzero];
csvwrite('ImpossibleDataset.csv', data);