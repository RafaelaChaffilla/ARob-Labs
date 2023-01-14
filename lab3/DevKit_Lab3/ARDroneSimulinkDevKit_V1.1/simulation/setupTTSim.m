%% =======================================================================
%  ARDrone Simulation Example: Trajectory Tracking 
%  ========================================================================
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
%bdclose all;
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
disp('    (3) Circle - Constant Yaw'); 
disp('    (4) Circle - Variable Yaw'); 
disp('    (5) Line of Sight'); 

choice = input('');

psi_ki = 0.1;
psi_kp = -1;
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
            plot(sim_results(c).tout, sim_results(c).h_sim.signals.values);
        end
        xlim([0 30]);
        legend('Pedida', sprintf('k = %0.1f', k_w_vector(1)), ...
            sprintf('k_w = %0.1f', k_w_vector(2)), ...
            sprintf('k_w = %0.1f', k_w_vector(3)), ...
            sprintf('k_w = %0.1f', k_w_vector(4)), ...
            sprintf('k_w = %0.1f', k_w_vector(5)), ...
            sprintf('k_w = %0.1f', k_w_vector(6)), ...
            'location', 'southeast');
        xlabel('Tempo [s]');
        ylabel('Altitude h [m]');
        saveas(fig1, '../imgs/Altitude_Sim', 'png');

    case 2
        % Horizontal Line Simulation
        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q1 = diag([2; 2; 2; 20; 20; 20])*2;
        Q2 = diag([1; 1; 1; 2.5; 2.5; 2.5]);
        R = diag([1 1 1]);
        
        k_w = 0.1;
        K = lqr(A,B,Q2,R);

        sim_line_results1 = sim('ARDroneTTSim_2019');

        % Plot motion
        fig2 = figure();
        hold on;
        plot(sim_line_results1.pd.signals.values(1, :), sim_line_results1.pd.signals.values(2, :));
        plot(sim_line_results1.x_sim.signals.values, sim_line_results1.y_sim.signals.values);
        xlim([0 5]);
        ylim([0 5]);
        legend('Pedida', 'Simulada', 'location', 'northwest');
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        saveas(fig2, '../imgs/Horizontal_line_Sim', 'png');
        
        fig2_3d = figure();
        plot3(sim_line_results1.x_sim.signals.values, sim_line_results1.y_sim.signals.values, sim_line_results1.h_sim.signals.values);
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        zlabel('Posição z [m]');
        saveas(fig2_3d, '../imgs/Horizontal_line_Sim_3d', 'png');
        
    case 3
        % Circle Simulation: const yaw angle
        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q1 = diag([2; 2; 2; 20; 20; 20])*2;
        Q2 = diag([1; 1; 1; 2.5; 2.5; 2.5]);
        R = diag([1 1 1]);

        k_w = 1;
        K = lqr(A,B,Q1,R);

        sim_circ_results1 = sim('ARDroneTTSim_2019');
        
        k_w = 0.1;
        K = lqr(A,B,Q2,R);

        sim_circ_results2 = sim('ARDroneTTSim_2019');

        % Plot motion
        fig3 = figure();
        hold on;
        plot(sim_circ_results1.pd.signals.values(1, :), sim_circ_results1.pd.signals.values(2, :));
        plot(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values);
        plot(sim_circ_results2.x_sim.signals.values, sim_circ_results2.y_sim.signals.values);
        xlim([-2 0.5]);
        ylim([-1  1]);
        legend('Pedida', 'Simulada: K_1 e k_w = 1', 'Simulada: K_2 e k_w = 0.1', 'location', 'northeast');
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        saveas(fig3, '../imgs/Circ_line_Sim_1', 'png');
        
        fig3_3d = figure();
        plot3(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values, sim_circ_results1.h_sim.signals.values);
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        zlabel('Posição z [m]');
        saveas(fig3_3d, '../imgs/Circ_line_Sim_1_3d', 'png');
        
    case 4
        % Circle Simulation: variable yaw angle
        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q1 = diag([2; 2; 2; 20; 20; 20])*2;
        Q2 = diag([1; 1; 1; 2.5; 2.5; 2.5]);
        R = diag([1 1 1]);

        k_w = 1;
        K = lqr(A,B,Q1,R);

        sim_circ_results1 = sim('ARDroneTTSim_2019');
        
        k_w = 0.1;
        K = lqr(A,B,Q2,R);

        sim_circ_results2 = sim('ARDroneTTSim_2019');

        % Plot motion
        fig4 = figure();
        hold on;
        plot(sim_circ_results1.pd.signals.values(1, :), sim_circ_results1.pd.signals.values(2, :));
        plot(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values);
        plot(sim_circ_results2.x_sim.signals.values, sim_circ_results2.y_sim.signals.values);
        xlim([-2 0.5]);
        ylim([-1 1]);
        legend('Pedida', 'Simulada: K_1 e k_w = 1', 'Simulada: K_2 e k_w = 0.1', 'location', 'northwest');
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        saveas(fig4, '../imgs/Circ_line_Sim_2', 'png');
        
        fig5 = figure();
        hold on;
        plot(sim_circ_results1.psi_d.time, sim_circ_results1.psi_d.signals.values);
        plot(sim_circ_results1.yaw_rad.time, sim_circ_results1.yaw_rad.signals.values); %simulada
        xlim([20 60]);
        legend('Pedida', 'Simulada', 'location', 'northwest');
        xlabel('Tempo [s]');
        ylabel('\psi [rad]');
        saveas(fig5, '../imgs/Circ_line_psi_Sim', 'png');
        
        fig4_3d = figure();
        plot3(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values, sim_circ_results1.h_sim.signals.values);
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        zlabel('Posição z [m]');
        saveas(fig4_3d, '../imgs/Circ_line_Sim_2_3d', 'png');
        
    case 5
        % Line of Sight
        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q1 = diag([2; 2; 2; 20; 20; 20])*2;
        Q2 = diag([1; 1; 1; 2.5; 2.5; 2.5]);
        R = diag([1 1 1]);

        k_w = 1;
        K = lqr(A,B,Q1,R);
        K = [zeros(3), K(:,4:6)];

        sim_circ_results1 = sim('ARDroneTTSim_2019');
        
        k_w = 0.1;
        K = lqr(A,B,Q2,R);
        K = [zeros(3), K(:,4:6)];
        
        sim_circ_results2 = sim('ARDroneTTSim_2019');

        % Plot motion
        fig6 = figure();
        hold on;
        plot(sim_circ_results1.pd.signals.values(1, :), ones(size(sim_circ_results1.pd.signals.values(1, :))));
        plot(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values);
        plot(sim_circ_results2.x_sim.signals.values, sim_circ_results2.y_sim.signals.values);
        % xlim([0 30]);
        legend('Pedida', 'Simulada: K_1 e k_w = 1', 'Simulada: K_2 e k_w = 0.1', 'location', 'southeast');
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        saveas(fig6, '../imgs/Los_Sim', 'png');
        
        fig7 = figure();
        hold on;
        plot(sim_circ_results1.psi_d.time, sim_circ_results1.psi_d.signals.values);
        plot(sim_circ_results1.yaw_rad.time, sim_circ_results1.yaw_rad.signals.values); %simulada
        % xlim([0 30]);
        legend('Pedida', 'Simulada', 'location', 'northeast');
        xlabel('Tempo [s]');
        ylabel('\psi [rad]');
        saveas(fig7, '../imgs/Los_psi_Sim', 'png');
        
        fig6_3d = figure();
        plot3(sim_circ_results1.x_sim.signals.values, sim_circ_results1.y_sim.signals.values, sim_circ_results1.h_sim.signals.values);
        xlabel('Posição x [m]');
        ylabel('Posição y [m]');
        zlabel('Posição z [m]');
        saveas(fig6_3d, '../imgs/Los_Sim_3d', 'png');
        
    case 6
        % Psi Simulation for different k_i and k_p
        k_w = 1;
        K=zeros(3,6);
        psi_ki_vector = [0.1:0.4:0.9];
        psi_kp_vector = -[0.2:0.4:1];
        for i = 1:size(psi_ki_vector,2)
            psi_ki = psi_ki_vector(i);
            for j = 1:size(psi_kp_vector,2)
                psi_kp = psi_kp_vector(j);
                sim_results(j+(i-1)*size(psi_kp_vector,2)) = sim('ARDroneTTSim_2019');
            end
        end

        % Plot altitude tracking for different k_w
        fig7 = figure();
        hold on;
        plot(sim_results(1).psi_d.time, sim_results(1).psi_d.signals.values);
        for c = 1:size(sim_results,2)
            plot(sim_results(c).tout, sim_results(c).yaw_rad.signals.values);
        end
        legend('pedida', 'k_i = 0.1; k_p = -0.2',...
            'k_i = 0.1; k_p = -0.6', 'k_i = 0.1; k_p = -1',...
            'k_i = 0.5; k_p = -0.2',...
            'k_i = 0.5; k_p = -0.6', 'k_i = 0.5; k_p = -1',...
            'k_i = 0.9; k_p = -0.2',...
            'k_i = 0.9; k_p = -0.6', 'k_i = 0.5; k_p = -1',...
            'location', 'southeast');
        xlabel('Tempo [s]');
        ylabel('\psi [rad]');
        saveas(fig7, '../imgs/PSI_Sim', 'png');
        
        
        
    otherwise
        %Loading Simulink model of ARDrone
        ARDroneTTSim_2019 ;
        
end