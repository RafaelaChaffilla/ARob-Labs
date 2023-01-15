%% LOS Study
%
% This script studies the overshoot and distance to settle given the
% velocity V and distance to look \Delta.
%
%% Script - Preparation for simulation (like TTSim)
close all;
clear all;
clc
% Adding ARDrone library path 
addpath './lib'; 
% Flight management system sample time. This is the sample time at which
% the control law is executed. 
FMS.Ts = 0.03;
% Time delay due to communication between drone and host computer
timeDelay = FMS.Ts*4; 
% Loading state space representation of vehicle dynamics
setupARModel;
% Simulation time
simDT = 0.005 ;
% Controller Gains (for simulation, experimental were different)
A = [zeros(3,3) eye(3);...
     zeros(3,3) zeros(3,3)];
B = [zeros(3,3); eye(3)];
Q1 = diag([2; 2; 2; 20; 20; 20])*2;
R = diag([1 1 1]);
k_w = 1;
K = lqr(A,B,Q1,R);
psi_ki = 0.1;
psi_kp = -1;

%% Gathering Data
Ratio_Array     = [0.5;1;2;5];
V_array         = linspace(0.1,0.5,20);
disp(['Starting Study']);
for c=1:length(V_array)
    V = V_array(c);
    for d=1:length(Ratio_Array)
        Delta = V*Ratio_Array(d);
        DATA = sim('LOS_Study_Sim');
        % Index of the first approach
        [~,result]      = max(logical(DATA.D_Approach(:)),[],1);
        D_APPROACH(c,d) = DATA.D_Approach(result);
        % Index of the first overshoot
        [~,result]      = max(logical(DATA.Overshoot(:)),[],1);
        OVERSHOOT(c,d)  = max([(DATA.Overshoot(result) - 1)*100,0]);
    end
    disp(['V = ' num2str(V) ' done']);
end
disp(['Ending Study']);

%% APPROACH
figure();
hold on
for c=1:length(Ratio_Array)
    plot(V_array, D_APPROACH(:,c));
    LEGEND{c} = ['\Delta/V = ' num2str(Ratio_Array(c))];
end
xlabel('V [ms^{-1}]');
ylabel('d_{approach} [m]');
legend(LEGEND, 'Location','eastoutside');
%% OVERSHOOT
figure();
hold on
for c=1:length(Ratio_Array)
    plot(V_array, OVERSHOOT(:,c));
    LEGEND{c} = ['\Delta/V = ' num2str(Ratio_Array(c))];
end
xlabel('V [ms^{-1}]');
ylabel('overshoot [%]');
legend(LEGEND, 'Location','eastoutside');