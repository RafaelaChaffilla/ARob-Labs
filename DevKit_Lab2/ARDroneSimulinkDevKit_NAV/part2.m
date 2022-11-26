clear all;
close all;

% Dacc  = load("./data/D_accelerometer.mat");
% Dgyr  = load("./data/D_gyro.mat");
% Dest1 = load("./data/D_estimated1.mat");
% Dpitch = load("./data/D_pitch.mat");
% DpitchA = load("./data/D_pitch_a.mat");
Dacc  = load("accelerometer.mat");
Dgyr  = load("gyro.mat");
Dest1 = load("estimated1.mat");
Dpitch = load("pitch.mat");
DpitchA = load("pitch_a.mat");
%% plots
plot_estimation(Dpitch, Dest1, DpitchA);

%% functions
function plot_estimation(pitch, est1, pitcha)
    figure();
    hold on;
    plot(pitch.ans(1,:), pitch.ans(2,:));
    plot(est1.ans(1,:), est1.ans(2,:));
%   plot(pitch.ans(1,:), mod(pitcha.ans(2,:),zeros(1,size(pitcha.ans,2))+360));
    legend('sensor','estimation_us', 'estimation_a');

end