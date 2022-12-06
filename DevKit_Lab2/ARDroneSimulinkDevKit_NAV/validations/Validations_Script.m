%% Filter Validations
%
%   Duarte Ferro Lopes  , n 95783
%   Jose Medeiro        , n 95811
%   Rafaela Chaffilla   , n 95840
%
%% Initialization
clear all;
close all;
addpath '../Kalman_Filters'
sampleTime  = 0.05;
%% Setup of filters
% 1st filter - no bias

Q           = 2*10^(-2);
R           = 1*10^(-2);
Kalman_1    = Setup_Kalman_1(Q, R, sampleTime);

disp(['1st Kalman Filter gains are L_1 = ' num2str(Kalman_1.L)])
% 2nd filter - bias

Q           = [8*10^(-3),1.5*10^(-3)];
R           = 1*10^(-2);
Kalman_2    = Setup_Kalman_2(Q, R, sampleTime);

disp(['2nd Kalman Filter gains are L_1 = ' num2str(Kalman_2.L(1))...
      '; L_2 = ' num2str(Kalman_2.L(2))]);

% P. Batista filter

Q           = [ 0.05, 0.05, 0.05,...
                0.01, 0.01, 0.01];
R           = [ 0.05, 0.05, 0.05];
Kalman_OP   = Setup_Kalman_OP(Q, R);

%% Tests
tests = 0;
%
% 1st filter
%   filter choise
filter = 1;
%   bias choise
bias = [0; 0; 0]*pi/180;
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%   bias choise
bias = [-1; 2; 0]*pi/180;
%   simulation
tests = tests +1;
SIM(2) = sim('Validation_Kalman_OP');

%
% 2nd filter
%   filter choise
filter = 2;
%   bias choise
bias = [5; -3; 0]*pi/180;
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_OP');

%
% 3nd filter
%   filter choise
filter = 3;
%   bias choise
bias = [5; -3; 3]*pi/180;
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_OP');

%% Plots