%% Plots Simulation

% xy plot
figure(1)
plot(sim_data.Data(1,:), sim_data.Data(2,:), sim_data.Data(7,:), sim_data.Data(8,:));
axis([0 5 -3 2]);
legend('desejado','actual')
xlabel('x [m]')
ylabel('y [m]')
print('traj_sim_exp4','-dpng');

% x plot
figure(2)
plot(sim_data.Time, sim_data.Data(1,:), sim_data.Time, sim_data.Data(7,:));
axis([0 25 0 5]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('x [m]')
print('x_sim_exp4','-dpng');

% y plot
figure(3)
plot(sim_data.Time, sim_data.Data(2,:), sim_data.Time, sim_data.Data(8,:));
axis([0 25 -2 2]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('y [m]')
print('y_sim_exp4','-dpng');

% z plot
figure(4)
plot(sim_data.Time, -sim_data.Data(3,:), sim_data.Time, sim_data.Data(9,:));
axis([0 25 0 2]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('z [m]')
print('z_sim_exp4','-dpng');

% phi plot
figure(5)
plot(sim_data.Time, 180/pi*sim_data.Data(4,:), sim_data.Time, 180/pi*sim_data.Data(10,:));
axis([0 25 -10 10]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('\phi [rad]')
print('phi_sim_exp4','-dpng');

% theta plot
figure(6)
plot(sim_data.Time, 180/pi*sim_data.Data(5,:), sim_data.Time, 180/pi*sim_data.Data(11,:));
axis([0 25 -10 10]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('\theta [rad]')
print('theta_sim_exp4','-dpng');

% psi plot
figure(7)
plot(sim_data.Time, 180/pi*wrapToPi(sim_data.Data(6,:)), sim_data.Time, 180/pi*wrapToPi(sim_data.Data(12,:)));
axis([0 25 -180 180]);
legend('desejado','actual')
xlabel('t [s]')
ylabel('\psi [rad]')
print('psi_sim_exp4','-dpng');


% p error plot
figure(8)
plot(sim_data.Time, sim_data.Data(13,:), sim_data.Time, sim_data.Data(14,:), sim_data.Time, sim_data.Data(15,:)-2);
axis([0 25 -1 1]);
legend('x','y','z')
xlabel('t [s]')
ylabel('erro [m]')
print('erro_posicoes_sim_exp4','-dpng');

% angle error plot
figure(9)
plot(sim_data.Time, 180/pi*sim_data.Data(16,:), sim_data.Time, 180/pi*sim_data.Data(17,:), sim_data.Time, 180/pi*wrapToPi(sim_data.Data(18,:)));
axis([0 25 -10 10]);
legend('\phi','\theta','\psi')
xlabel('t [s]')
ylabel('erro [º]')
print('erro_angulos_sim_exp4','-dpng');



%% Plots Real

% figure(1)
% plot(states.time(1:3065), 180/pi*states.signals.values(802:3866,1));
% axis([0 31 -10 10]);
% xlabel('t [s]')
% ylabel('\phi [º]')
% print('phi_exp3','-dpng');
% 
% figure(2)
% plot(states.time(1:3065), 180/pi*states.signals.values(802:3866,2));
% axis([0 31 -10 10]);
% xlabel('t [s]')
% ylabel('\theta [º]')
% print('theta_exp3','-dpng');
% 
% figure(3)
% plot(states.time(1:3065), 180/pi*states.signals.values(802:3866,3));
% axis([0 31 -180 180]);
% xlabel('t [s]')
% ylabel('\psi [º]')
% print('psi_exp3','-dpng');
% 
% figure(4)
% plot(states.time(1:3065), states.signals.values(802:3866,4));
% axis([0 31 -1 1]);
% xlabel('t [s]')
% ylabel('u [m/s]')
% print('u_exp3','-dpng');
% 
% figure(5)
% plot(states.time(1:3065), states.signals.values(802:3866,5));
% axis([0 31 -1 1]);
% xlabel('t [s]')
% ylabel('v [m/s]')
% print('v_exp3','-dpng');
% 
% figure(6)
% plot(states.time(1:3065), states.signals.values(802:3866,6));
% axis([0 31 -1 1]);
% xlabel('t [s]')
% ylabel('w [m/s]')
% print('w_exp3','-dpng');
% 
% figure(7)
% plot(states.time(1:3065), states.signals.values(802:3866,7));
% axis([0 31 -3 1]);
% xlabel('t [s]')
% ylabel('x [m]')
% print('x_exp3','-dpng');
% 
% figure(8)
% plot(states.time(1:3065), states.signals.values(802:3866,8));
% axis([0 31 -2 2]);
% xlabel('t [s]')
% ylabel('y [m]')
% print('y_exp3','-dpng');
% 
% figure(9)
% plot(states.time(1:3065), states.signals.values(802:3866,9));
% axis([0 31 0 2]);
% xlabel('t [s]')
% ylabel('z [m]')
% print('z_exp3','-dpng');
% 
% figure(10)
% plot(states.signals.values(802:3866,7), states.signals.values(802:3866,8));
% axis([-3 1 -2 2]);
% xlabel('x [m]')
% ylabel('y [m]')
% print('traj_exp3','-dpng');