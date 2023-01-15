%% PLOT_EXP
%
% This script plots the experimental data obtained during the tests and
% present in /data
%
%%
close all;
clear all;

PATH_DAT = './data/';
PATH_IMG = './img/';
% Results : 
%   Line 1   : Time
%   Line 2-10: States
%
% Desired : 
%   Line 1   : Time
%   Line 2-4 : P_d
%   Line 5   : Psi_d
%   Line 6   : Roll_ref
%   Line 7   : Pitch_ref
%   Line 8   : Commands enable Flag

%% Altitude
Alt_exp.desired = load([PATH_DAT 'desired_data_1.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_1.mat']);

FIG = figure();
hold on;
plot(Alt_exp.desired.ans(1,:), -Alt_exp.desired.ans(4,:));
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(10,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Tempo [s]');
ylabel('Altitude h [m]');
saveas(FIG,[PATH_IMG 'Altitude_EXP'], 'png');

%% Trajétória Horizontal - V = 0.2m/s, 45deg
Alt_exp.desired = load([PATH_DAT 'desired_data_2.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_2.mat'    ]);

FIG = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northwest');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
saveas(FIG, [PATH_IMG 'Horz_45_EXP'], 'png');

FIG = figure();
subplot(1,3,1);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(7,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\theta [º]');
xlim([21 36]);

subplot(1,3,2);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(6,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(2,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\phi [º]');
xlim([21 36]);

subplot(1,3,3);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(5,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\psi [º]');
xlim([21 36]);
saveas(FIG, [PATH_IMG 'Horz_45_EXP_pitch'], 'png');

%% Trajétória Horizontal - V = 0.2m/s, 0deg
Alt_exp.desired = load([PATH_DAT 'desired_data_2_2.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_2_2.mat']);

FIG = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:) );
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:) );
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-0.4 0.4]);
saveas(FIG, [PATH_IMG 'Horz_0_EXP'], 'png');

FIG = figure();
subplot(1,3,1);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(7,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\theta [º]');
xlim([10 80]);

subplot(1,3,2);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(6,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(2,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\phi [º]');
xlim([10 80]);

subplot(1,3,3);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(5,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\psi [º]');
xlim([10 80]);
saveas(FIG, [PATH_IMG 'Horz_0_EXP_pitch'], 'png');

%% Trajétória Circular - Psi Constante
Alt_exp.desired = load([PATH_DAT 'desired_data_3.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_3.mat']);

FIG = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-1.2 1.2]);
xlim([-2.2 0.7]);
saveas(FIG, [PATH_IMG 'circ_yawcte_EXP'], 'png');

FIG = figure();
subplot(1,3,1);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(7,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\theta [º]');
xlim([10 80]);

subplot(1,3,2);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(6,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(2,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\phi [º]');
xlim([10 80]);

subplot(1,3,3);
hold on;
plot(Alt_exp.desired.ans(1,:), Alt_exp.desired.ans(5,:)*180/pi);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(3,:)*180/pi);
legend('Pedida', 'Experimental', 'location', 'southwest');
xlabel('Tempo [s]');
ylabel('\psi [º]');
xlim([10 80]);

saveas(FIG, [PATH_IMG 'circ_yawcte_EXP_pitch'], 'png');
%% Trajétória Circular - Psi Variável
Alt_exp.desired = load([PATH_DAT 'desired_data_4.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_4.mat']);

% WHEN COMMANDS WERE ENABLED
[~,Enabled] 	= max(logical(Alt_exp.desired.ans(8,:)'),[],1);
[~,Active]      = max(logical(Alt_exp.desired.ans(5,:)'),[],1);

% ANGLES AFTER ACTIVE
TIME    = Alt_exp.desired.ans(1,Active:end);
PSI.DES = Alt_exp.desired.ans(5,Active:end);
PSI.EXP = Alt_exp.results.ans(4,Active:end);
PHI.DES = Alt_exp.desired.ans(6,Active:end);
PHI.EXP = Alt_exp.results.ans(2,Active:end);
THE.DES = Alt_exp.desired.ans(7,Active:end);
THE.EXP = Alt_exp.results.ans(3,Active:end);
% X AND Y
X.DES   = Alt_exp.desired.ans(2,Active:end);
X.EXP   = Alt_exp.results.ans(8,Active:end);
Y.DES   = Alt_exp.desired.ans(3,Active:end);
Y.EXP   = Alt_exp.results.ans(9,Active:end);
% DEMODDING OF PSI
counter = 0;
for i=2:length(PSI.EXP)
    if PSI.EXP(i)+counter - PSI.EXP(i-1) > pi
        counter = counter - 2*pi;
    end
    if PSI.EXP(i)+counter - PSI.EXP(i-1) < -pi
        counter = counter + 2*pi;
    end
    PSI.EXP(i) = PSI.EXP(i) + counter;
end
counter = 0;
for i=2:length(PSI.DES)
    if PSI.DES(i)+counter - PSI.DES(i-1) > pi
        counter = counter - 2*pi;
    end
    if PSI.DES(i)+counter - PSI.DES(i-1) < -pi
        counter = counter + 2*pi;
    end
    PSI.DES(i) = PSI.DES(i) + counter;
end

% PLOT POSISAO
FIG = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-1 1]);
xlim([-2 0.5]);
saveas(FIG, [PATH_IMG 'YAWVAR_PATH_EXP'], 'png');

% ERRORS
FIG = figure();
hold on;
plot(TIME, X.DES);
plot(TIME, X.EXP);
ylabel('Posição x [m]');
xlabel('Tempo [s]');
legend('Pedido', 'Experimental', 'location', 'southeast');
saveas(FIG, [PATH_IMG 'YAWVAR_X_EXP'], 'png');
FIG = figure();
hold on;
plot(TIME, Y.DES);
plot(TIME, Y.EXP);
legend('y');
ylabel('Posição y [m]');
xlabel('Tempo [s]');
legend('Pedido', 'Experimental', 'location', 'southeast');
saveas(FIG, [PATH_IMG 'YAWVAR_T_EXP'], 'png');

% PLOT ANGLES
FIG = figure();
subplot(1,3,1);
hold on;
plot(TIME,(PSI.EXP - PSI.DES)*180/pi);
xlim([TIME(Active) 50]);
legend('Erro \psi', 'location', 'southeast');
xlabel('Tempo [s]');
ylabel('e [º]');
subplot(1,3,2);
hold on;
plot(TIME,THE.DES*180/pi);
plot(TIME,THE.EXP*180/pi);
xlim([TIME(Active) 50]);
legend('Pedida', 'Experimental', 'location', 'southeast');
xlabel('Tempo [s]');
ylabel('\theta [º]');
subplot(1,3,3);
hold on;
plot(TIME, PHI.DES*180/pi);
plot(TIME, PHI.EXP*180/pi);
xlim([TIME(Active) 50]);
legend('Pedida', 'Experimental', 'location', 'southeast');
xlabel('Tempo [s]');
ylabel('\phi [º]');
saveas(FIG, [PATH_IMG 'YAWVAR_ANG_EXP'], 'png');

%% LOS - first try
Alt_exp.desired = load([PATH_DAT 'desired_data_5.mat']);
Alt_exp.results = load([PATH_DAT 'exp_data_5.mat']);

% WHEN COMMANDS WERE ENABLED
[~,Enabled] 	= max(logical(Alt_exp.desired.ans(8,:)'),[],1);
[~,Active]      = max(logical(Alt_exp.desired.ans(5,:)'),[],1);
% DISTANCE SINCE ENABLED
x_0	= Alt_exp.results.ans(8,Enabled);
y_0	= Alt_exp.results.ans(9,Enabled);

% PSI AFTER ENABLED
QUIVER_D = downsample(Alt_exp.desired.ans(:,Active:end)',30);
QUIVER_D = QUIVER_D';
QUIVER_R = downsample(Alt_exp.results.ans(:,Active:end)',30);
QUIVER_R = QUIVER_R';

FIG = figure();
hold on;
plot([-0.4, 1.4],[y_0 + 1, y_0 + 1]);
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
plot(x_0,y_0,'*r');
quiver(QUIVER_R(8,:), QUIVER_R(9,:),...
       cos(QUIVER_R(4,:)), sin(QUIVER_R(4,:)),0.5);
legend('Pedida', 'Experimental',...
       'Inicial','\psi',...
       'location', 'southeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-0.2 1.4]);
xlim([-0.4 1]);
saveas(FIG, [PATH_IMG 'LOS_PATH_EXP'], 'png');
% ANGLES
psi_d = mod(Alt_exp.desired.ans(5,:), 2*pi) ;
for i=1:size(psi_d,2)
    if psi_d(i) > pi
        psi_d(i) = -2*pi + psi_d(i) ; 
    end
end

FIG = figure();
hold on;
plot(Alt_exp.desired.ans(1,:), psi_d);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(4,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Tempo [s]');
ylabel('\Psi [rad]');
xlim([9.12 32]);
saveas(FIG, [PATH_IMG 'LOS_YAW_EXP'], 'png');