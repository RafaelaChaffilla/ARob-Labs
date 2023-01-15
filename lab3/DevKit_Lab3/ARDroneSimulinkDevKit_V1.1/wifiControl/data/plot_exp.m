close all;
clear all;

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
Alt_exp.desired = load('desired_data_1.mat');
Alt_exp.results = load('exp_data_1.mat');

fig1 = figure();
hold on;
plot(Alt_exp.desired.ans(1,:), -Alt_exp.desired.ans(4,:));
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(10,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Tempo [s]');
ylabel('Altitude h [m]');
saveas(fig1, '../../imgs/Altitude_EXP', 'png');

%% Trajétória Horizontal - V = 0.2m/s, 45deg
Alt_exp.desired = load('desired_data_2.mat');
Alt_exp.results = load('exp_data_2.mat');

fig2 = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northwest');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
saveas(fig2, '../../imgs/Horz_45_EXP', 'png');

fig2a = figure();
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
saveas(fig2a, '../../imgs/Horz_45_EXP_pitch', 'png');

%% Trajétória Horizontal - V = 0.2m/s, 0deg
Alt_exp.desired = load('desired_data_2_2.mat');
Alt_exp.results = load('exp_data_2_2.mat');

fig3 = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:) );
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:) );
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-0.4 0.4]);
saveas(fig3, '../../imgs/Horz_0_EXP', 'png');

fig3a = figure();
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
saveas(fig3a, '../../imgs/Horz_0_EXP_pitch', 'png');

%% Trajétória Circular - Psi Constante
Alt_exp.desired = load('desired_data_3.mat');
Alt_exp.results = load('exp_data_3.mat');

fig4 = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-1.2 1.2]);
xlim([-2.2 0.7]);
saveas(fig4, '../../imgs/circ_yawcte_EXP', 'png');

fig4a = figure();
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

saveas(fig4a, '../../imgs/circ_yawcte_EXP_pitch', 'png');
%% Trajétória Circular - Psi Variável
Alt_exp.desired = load('desired_data_4.mat');
Alt_exp.results = load('exp_data_4.mat');

fig5 = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), Alt_exp.desired.ans(3,:));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
ylim([-1 1]);
xlim([-2 0.5]);
saveas(fig5, '../../imgs/circ_yawvar_EXP', 'png');

psi_d = mod(Alt_exp.desired.ans(5,:), 2*pi) ;
for i=1:size(psi_d,2)
    if psi_d(i) > pi
        psi_d(i) = -2*pi + psi_d(i) ; 
    end
end

fig8 = figure();
hold on;
plot(Alt_exp.desired.ans(1,:), psi_d);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(4,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Tempo [s]');
ylabel('\Psi [rad]');
xlim([17.22 55]);
saveas(fig8, '../../imgs/yawvar_EXP', 'png');

%% LOS - first try
Alt_exp.desired = load('desired_data_5.mat');
Alt_exp.results = load('exp_data_5.mat');

fig6 = figure();
hold on;
plot(Alt_exp.desired.ans(2,:), ones(size(Alt_exp.desired.ans(2,:))));
plot(Alt_exp.results.ans(8,:), Alt_exp.results.ans(9,:));
legend('Pedida', 'Experimental', 'location', 'southeast');
xlabel('Posição x [m]');
ylabel('Posição y [m]');
% ylim([-1 1]);
% xlim([-2 0.5]);
saveas(fig6, '../../imgs/LOS1_EXP', 'png');

psi_d = mod(Alt_exp.desired.ans(5,:), 2*pi) ;
for i=1:size(psi_d,2)
    if psi_d(i) > pi
        psi_d(i) = -2*pi + psi_d(i) ; 
    end
end

fig7 = figure();
hold on;
plot(Alt_exp.desired.ans(1,:), psi_d);
plot(Alt_exp.results.ans(1,:), Alt_exp.results.ans(4,:));
legend('Pedida', 'Experimental', 'location', 'northeast');
xlabel('Tempo [s]');
ylabel('\Psi [rad]');
xlim([9.12 32]);
saveas(fig7, '../../imgs/los_yaw_EXP', 'png');