%% =======================================================================
%  Start Here Script

%  =======================================================================
%  bench
%  This script guides the user through examples of both simulations and
%  Wi-Fi control models for the Parrot ARDrone.
%  Authors:
%       David Escobar Sanabria -> descobar@aem.umn.edu
%       Pieter J. Mosterman -> pieter.mosterman@mathworks.com
%  =======================================================================

%%
%  Cleaning workspace
bdclose all;
clear all;
bdclose all;
clc
sampleTime=0.05;
%% Kalman Filters Setup
addpath './Kalman_Filters/';
Q           = 1*10^(-2);
R           = 1*10^(-2);
Kalman_1    = Setup_Kalman_1(Q,R,sampleTime);

Q           = [6*10^(-2),1.5*10^(-2)];
R           = 5*10^(-2);
Kalman_2    = Setup_Kalman_2(Q,R,sampleTime);

Q           = [ ones(3,1)*4*10^(-3);...
                ones(3,1)*1*10^(-8);];
R           =   ones(3,1)*8*10^(-3);
Kalman_OP   = Setup_Kalman_OP(Q,R);
%% Selection

disp('This script guides the user through examples of both simulations and');
disp('Wi-Fi control Simulink models for the Parrot ARDrone ');
disp(' '); 

disp('Select one of the following options:'); 
disp('    (1) Simulation'); 
disp('    (2) Wi-Fi control. The computer shoud be connected to the drone'); 
disp('    (3) Replay from stored data'); 

option = input('');


%%
switch option
    case 1
        
        disp('Select one of the following options for simulation:'); 
        disp('    (1) Baseline simulation: The ARDrone block with inputs and scopes to visualize outputs'); 
        disp('    (2) Hover: Vehicle is held at constant position'); 
        disp('    (3) Waypoint tracking: Vehicle tracks a list of waypoints'); 
        option2 = input(' ');
        
        switch option2
            case 1
                cd simulation; 
                setupBaseModel; 
            case 2
                cd simulation; 
                setupHoverSim; 
            case 3
                cd simulation; 
                setupWPTrackingSim;
            otherwise
                disp('An incorrect option was selected')
                
        end
        
        
    case 2
        disp('Select one of the following options for Wi-Fi control:'); 
        disp('    (1) Hover: Vehicle is held at constant position'); 
        disp('    (2) Waypoint tracking: Vehicle tracks a list of waypoints'); 
        option2 = input(' ');
        
        switch option2
            case 1
                cd wifiControl; 
                setupHover; 
                % Building model using RTWT. Install rtwt if not installed 
                % using 'rtwintgt -setup'
                % rtwbuild('ARDroneHover_V2');
            case 2
                cd wifiControl; 
                setupWPTracking;  
                % Building model using RTWT. Install rtwt if not installed 
                % using 'rtwintgt -setup'
                % rtwbuild('ARDroneWPTracking_V2');
            otherwise
                disp('An incorrect option was selected')
        end
    case 3
        disp('Filename from where navdata will be loaded:'); 
        filename = input(' ','s');
        loadedNavdata = load(filename);
        if isfield(loadedNavdata,'ans')
            loadedNavdata = loadedNavdata.ans;
            finalTime = loadedNavdata.Time(end);
        else
            loadedNavdata = loadedNavdata.navdata;
            finalTime = loadedNavdata.time(end);
        end
        cd Replay;
        setupReplay;
    otherwise
       disp('An incorrect option was selected')

end





