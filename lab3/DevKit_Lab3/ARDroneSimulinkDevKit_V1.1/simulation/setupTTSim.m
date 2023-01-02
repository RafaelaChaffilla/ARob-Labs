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
%   Modified by Group 2:
%
%   Duarte Lopes,      95783
%   José Medeiro,      95811
%   Rafaela Chaffilla, 95840
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

%% choose simulation to run

disp('Choose the simulation you want to run:'); 
disp('    (1) Altitude'); 
disp('    (2) Horizontal Line'); 

choice = input('');

switch choice
    case 1
        % Altitude Simulation for different k_w
        k_w_vector=[0.5:0.5:3];
        K=zeros(3,6);
        for i = 1:size(k_w_vector,2)
            k_w = k_w_vector(i);
            sim_results(i) = sim('ARDroneTTSim_2019');
        end

        % Plot altitude tracking for different k_w
        fig1 = figure();
        hold on;
        plot(sim_results(1).pd.time, sim_results(1).pd.signals.values(3, :));
        for c = 1:size(sim_results,2)
            plot(sim_results(c).tout, - sim_results(c).h_sim.signals.values);
        end
        xlim([0 30]);
        legend('pedida', sprintf('k = %0.1f', k_w_vector(1)), ...
            sprintf('k_w = %0.1f', k_w_vector(2)), ...
            sprintf('k_w = %0.1f', k_w_vector(3)), ...
            sprintf('k_w = %0.1f', k_w_vector(4)), ...
            sprintf('k_w = %0.1f', k_w_vector(5)), ...
            sprintf('k_w = %0.1f', k_w_vector(6)), ...
            'location', 'west');
        xlabel('Tempo [s]');
        ylabel('Altitude h [m]');
        saveas(fig1, '../imgs/Altitude_Sim', 'png');

    case 2
        % Horizontal Line Simulation
        k_w = 1;

        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q = diag([2; 2; 2; 20; 20; 20]);
        R = diag([1 1 1]);

        K = lqr(A,B,Q,R);

        sim_line_results = sim('ARDroneTTSim_2019');

        % Plot altitude tracking for different k_w
        fig2 = figure();
        hold on;
        plot(sim_line_results.pd.signals.values(1, :), sim_line_results.pd.signals.values(2, :));
        plot(sim_line_results.x_sim.signals.values, sim_line_results.y_sim.signals.values);
        % xlim([0 30]);
        legend('pedida', 'simulada', 'location', 'northwest');
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        saveas(fig2, '../imgs/Horizontal_line_Sim', 'png');
        
    otherwise
        %Loading Simulink model of ARDrone
        ARDroneTTSim ;
        
end


 
 
 
 
