clear all;
close all;

Dacc  = load("./data/D_accelerometer.mat");
Dest1 = load("./data/D_estimated1.mat");
% Dpitch = load("./data/D_pitch.mat");
% DpitchA = load("./data/D_pitch_a.mat");
%Dacc  = load("accelerometer.mat");
%Dgyr  = load("gyro.mat");
Dest1 = load("estimated1.mat");
Dpitch = load("pitch.mat");
%% plots
plot_estimation(Dpitch, Dest1);

%% functions
function plot_estimation(pitch, est1)
    figure();
    hold on;
    plot(pitch.ans(1,:), pitch.ans(2,:));
    plot(est1.ans(1,:), est1.ans(2,:));
    legend('sensor','estimation');

end