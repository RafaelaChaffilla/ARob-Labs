close all;
clear all;

% Results : 
%   Line 1   : Time
%   Line 2-10: States
%
% Desired :
%   Line 1-3 : P_d
%   Line 4   : Psi_d
%   Line 5   : Roll_ref
%   Line 6   : Pitch_ref
%   Line 7   : Commands enable Flag

%% Altitude_Old
Alt_exp.desired = load('desired_data_1.mat');
Alt_exp.results = load('exp_data_1.mat');

fig1 = figure();
hold on;
plot(Alt_exp.desired.ans.time, Alt_exp.desired.ans.data(3, :));
plot(Alt_exp.results.ans.h__m_.time, Alt_exp.results.ans.h__m_.data);
legend('Pedida', 'Experimental', 'location', 'northwest');
xlabel('Tempo [s]');
ylabel('Altitude h [m]');
xlim([26 48]);
saveas(fig1, '../../../imgs/Altitude_EXP_OG', 'png');