%% =======================================================================
%  ARDrone Simulation Example: Trajectory Tracking 
%  =======================================================================
%  
%  The simulation is used to validate the trajectory tracking controller
%  and generation of the desired trajectory for the ARDRone before flight 
%  testing. Control and guidance blocks are
%  exactly the same for both simulation and real-time Wi-Fi control.
%  
%  =======================================================================

%%
%  Cleaning workspace
bdclose all;
close all;
clear all;
clc

%%
% Adding ARDrone library path 
addpath ../lib; 
%% Simulation parameters

% Flight management system sample time. This is the sample time at which
% the control law is executed. 
FMS.Ts = 0.03;%0.065; 

% Time delay due to communication between drone and host computer
timeDelay = FMS.Ts*4; 


%% Vehicle model based on linear dynamics

% Loading state space representation of vehicle dynamics
setupARModel; 

%%
% Loading list of waypoints
% waypoints = getWaypoints() ;


%% 
% Simulation time
simDT = 0.005 ;

%% Altitude Simulation for different k_w
% k_w_vector=[0.5:0.5:3];
% for i = 1:size(k_w_vector,2)
%     k_w = k_w_vector(i);
%     sim_results(i) = sim('ARDroneTTSim');
% end
% 
% % Plot altitude tracking for different k_w
% figure();
% hold on;
% plot(sim_results(1).pd.time, sim_results(1).pd.signals.values(3, :));
% for c = 1:size(sim_results,2)
%     plot(sim_results(c).tout, - sim_results(c).h_sim.signals.values);
% end
% xlim([0 30]);
% legend('desired', sprintf('k = %0.1f', k_w_vector(1)), ...
%     sprintf('k = %0.1f', k_w_vector(2)), ...
%     sprintf('k = %0.1f', k_w_vector(3)), ...
%     sprintf('k = %0.1f', k_w_vector(4)), ...
%     sprintf('k = %0.1f', k_w_vector(5)), ...
%     sprintf('k = %0.1f', k_w_vector(6)), ...
%     'location', 'west');
% xlabel('Tempo [s]');
% ylabel('Altitude h [m]');

%% Horizontal Line Simulation
k_w = 1.5;

A = [zeros(3,3) eye(3);...
     zeros(3,3) zeros(3,3)];
B = [zeros(3,3); eye(3)];

Q = diag([2; 2; 2; 20; 20; 20]);
R = diag([1 1 1]);

K = lqr(A,B,Q,R);

sim_line_results = sim('ARDroneTTSim');

% Plot altitude tracking for different k_w
figure();
hold on;
plot(sim_line_results.pd.signals.values(1, :), sim_line_results.pd.signals.values(2, :));
plot(sim_line_results.x_sim.signals.values, sim_line_results.y_sim.signals.values);
% xlim([0 30]);
legend('desired', 'simulated', 'location', 'west');
xlabel('Posi��o x [m]');
ylabel('Posi��o y [m]');

% Loading Simulink model of ARDrone
% ARDroneTTSim ;
 
 
 
 
