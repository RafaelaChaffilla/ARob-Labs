%% ========================================================================
%  ARDrone Wi-Fi Control: Hovering and Position Control
%  ========================================================================
%
%  This script initializes a Simulink model that enables control of the
%  Parrot ARDrone via Wi-Fi using Simulink. The control law 
%  tracks a  desired target position. 
%
%  IMPORTANT NOTE: The position estimation given in this example is inaccurate and
%  may lead to unexpected behavior of the vehicle.  The position is estimated
%  by intergating a velocity estimation, which in turn is inaccurate.
%  The velocity is estimated by the Drone onboard flight computer. The
%  velocity estimation can be improved if there are features on the floor
%  as the Drone runs an optical flow algorithm for velocity determination.
%  
%  
%  Requirements:
%       Matlab, Simulink, Real-Time-Wndows-Target
%
%  Authors:
%       David Escobar Sanabria -> descobar@aem.umn.edu
%       Pieter J. Mosterman -> pieter.mosterman@mathworks.com
%          
%  ========================================================================


%% Cleaning the workspace

bdclose all;
clear all;
close all; 

%% Adding ARDrone library to the path
addpath ../lib/ ; 

%%
%  Sample time of Simulink model. 
sampleTime = 0.03;%0.065; 

%%
%% choose simulation to run

disp('Choose the simulation you want to run:'); 
disp('    (1) Altitude'); 
disp('    (2) Horizontal Line'); 
disp('    (3) Circle - Constant Yaw'); 
disp('    (4) Circle - Variable Yaw'); 

choice = input('');

%% choose simulation to run

disp('Choose the simulation you want to run:'); 
disp('    (1) Altitude'); 
disp('    (2) Horizontal Line'); 
disp('    (3) Circle - Constant Yaw'); 
disp('    (4) Circle - Variable Yaw'); 
disp('    (5) LOS');

choice = input('');

switch choice
    case 1
        % Altitude Simulation for different k_w
        k_w=1;
        K=zeros(3,6);

    case 2
        % Horizontal Line Simulation
        k_w = 1;

        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q = diag([2; 2; 2; 20; 20; 20]);
        R = diag([1 1 1]);

        K = lqr(A,B,Q,R);
        
    case 3
        % Circle Simulation: const yaw angle
        k_w = 1;

        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q = diag([2; 2; 2; 20; 20; 20]).*2;
        R = diag([1 1 1]);

        K = lqr(A,B,Q,R);
%         K=[diag([4 4 4]) diag([6 6 6])];
        
    case 4
        % Circle Simulation: variable yaw angle
        k_w = 1;

        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q = diag([2; 2; 2; 20; 20; 20]).*2;
        R = diag([1 1 1]);

        K = lqr(A,B,Q,R);
%         K=[diag([4 4 4]) diag([6 6 6])];
     case 5
        % LOS
        k_w = 1;

        A = [zeros(3,3) eye(3);...
             zeros(3,3) zeros(3,3)];
        B = [zeros(3,3); eye(3)];

        Q = diag([2; 2; 2; 20; 20; 20]).*2;
        R = diag([1 1 1]);

        K = lqr(A,B,Q,R);
%         K=[diag([4 4 4]) diag([6 6 6])];
        
    otherwise
        %Loading Simulink model of ARDrone
        quit;
        
end

%ARDroneTT ; 