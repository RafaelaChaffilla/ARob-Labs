% Attention:
% This code only works in sections, so to save file X, one must run only
% section X

% Times when the experiment was being held

% A - from 0 to 30 s
% B - from 7.75 to 37.75 s
% C - from 10.035 to 40.035 s

%% A
clear all;
sampleTime = 0.05;
loadd = load('./data/A_navdata.mat');
loadedNavdata = loadd.navdata;
cd Replay;
finalTime  = 30;

set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchFile','FileName',sprintf('../dataA/A_pitch.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchAFile','FileName',sprintf('../dataA/A_pitch_a.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/accFile','FileName',sprintf('../dataA/A_acc.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/gyroFile','FileName',sprintf('../dataA/A_gyro.mat'));
setupReplay;

%% B
clear all;
sampleTime = 0.05;
loadd = load('./data/B_navdata.mat');
loadedNavdata = loadd.navdata;
cd Replay;
finalTime  = 37.75;

set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchFile','FileName',sprintf('../dataB/B_pitch.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchAFile','FileName',sprintf('../dataB/B_pitch_a.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/accFile','FileName',sprintf('../dataB/B_acc.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/gyroFile','FileName',sprintf('../dataB/B_gyro.mat'));
setupReplay;

%% C
clear all;
sampleTime = 0.05;
loadd = load('./data/C_navdata.mat');
loadedNavdata = loadd.navdata;
cd Replay;
finalTime  = 40.035;

set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchFile','FileName',sprintf('../dataC/C_pitch.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/pitchAFile','FileName',sprintf('../dataC/C_pitch_a.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/accFile','FileName',sprintf('../dataC/C_acc.mat'));
set_param('ARDroneReplay_V2/ARDrone Wi-Fi  Block NAV Replay/gyroFile','FileName',sprintf('../dataC/C_gyro.mat'));
setupReplay;

%% Remove first values
%A already starts at 0 seconds
%% B
clear all;
B_acc = load('./dataB/B_acc.mat');
B_gyr = load('./dataB/B_gyro.mat');
% find the point
t0 = find(B_acc.ans(1,:) == 7.75)
% selection
B_acc.ans(:,1:t0) = [];
B_gyr.ans(:,1:t0) = [];

save('./dataB/B_acc.mat','B_acc');
save('./dataB/B_gyro.mat','B_gyr');

%% C
clear all;
C_acc = load('./dataC/C_acc.mat');
C_gyr = load('./dataC/C_gyro.mat');
% find the point
%t0c = find(C_acc.ans(1,:) == 10.035);
t0c = 202;
% % selection
C_acc.ans(:,1:t0c) = [];
C_gyr.ans(:,1:t0c) = [];

% X=C_acc.C_acc;
% Y=C_gyr.C_gyr;
save('./dataC/C_acc.mat','C_acc');
save('./dataC/C_gyro.mat','C_gyr');

