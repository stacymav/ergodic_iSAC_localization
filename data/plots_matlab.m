clear all
close all

global param 
global master

%Read and format states and time
states = dlmread('states.csv');
time = states(:,1);
states(:,1) = [];


% Load all settings in struct "param" 
param = settings();


figure(1), plot(states(:,1))
figure(2), plot(states(:,2))
figure(3), plot(states(:,3))


figure(4), plot(states(:,4))
figure(5), plot(states(:,5))
figure(6), plot(states(:,6))

figure(7), plot(states(:,7), states(:,8))
axis equal
axis([0,1,0,1])

figure(2), plot(time,states(:,3))

%Plots
%plot(time,states(:,1:6))
%legend('x','y','z','phi','theta','psi')

% plot_control();