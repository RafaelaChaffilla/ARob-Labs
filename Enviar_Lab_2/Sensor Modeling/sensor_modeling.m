%% Filter Validations
%
%   Duarte Ferro Lopes  , n 95783
%   Jose Medeiro        , n 95811
%   Rafaela Chaffilla   , n 95840
%
%   This script uses the data that was treated previosly of experiences A,
%   B and C to obtain the sensors mean and variance values, as well as
%   plotting the acceleration and angular velocity obtained directly from
%   the sensors.
%% Initialization
clear all;
close all;

%% Loading of files
Aacc=load("./dataA/A_acc.mat");
Agyr=load("./dataA/A_gyro.mat");
Bacc=load("./dataB/B_acc.mat");
Bgyr=load("./dataB/B_gyro.mat");
Cacc=load("./dataC/C_acc.mat");
Cgyr=load("./dataC/C_gyro.mat");

%% plots
plot_acc_gyr(Aacc, Agyr);
plot_acc_gyr(Bacc.B_acc, Bgyr.B_gyr);
plot_acc_gyr(Cacc.C_acc, Cgyr.C_gyr);

%% Mean and Variance
% Choose valid data
Aacc.ans = Aacc.ans(:, 4:end);
Agyr.ans = Agyr.ans(:, 4:end);
Cacc.C_acc.ans = Cacc.C_acc.ans(:, 190:end);
Cgyr.C_gyr.ans = Cgyr.C_gyr.ans(:, 190:end);

[mean_Aacc, var_Aacc, mean_Agyr, var_Agyr] = mean_and_var(Aacc, Agyr);
[mean_Bacc, var_Bacc, mean_Bgyr, var_Bgyr] = mean_and_var(Bacc.B_acc, Bgyr.B_gyr);
[mean_Cacc, var_Cacc, mean_Cgyr, var_Cgyr] = mean_and_var(Cacc.C_acc, Cgyr.C_gyr);

%% Functions
function plot_acc_gyr(Aacc, Agyr)
    figure();
    hold on;
    for i=2:size(Aacc.ans,1)
        plot(Aacc.ans(1,:), Aacc.ans(i,:).*0.01); %units in m/s
    end
    xlim([Aacc.ans(1,1) Aacc.ans(1,end)]);
    legend('a_x','a_y','a_z');
    xlabel("Tempo [s]");
    ylabel("Acelaração Linear [m/s]");

    figure();
    hold on;
    for i=2:size(Agyr.ans,1)
        plot(Agyr.ans(1,:), Agyr.ans(i,:));
    end
    xlim([Aacc.ans(1,1) Aacc.ans(1,end)]);
    legend('p','q','r');
    xlabel("Tempo [s]");
    ylabel("Acelaração Angular [º/s]");

end

function [mean_acc, var_acc, mean_gyr, var_gyr] = mean_and_var(Aacc, Agyr)
    mean_acc = mean(Aacc.ans(2:4,:)')*0.01; %units em m/s
    mean_gyr = mean(Agyr.ans(2:4,:)');
    var_acc = var(Aacc.ans(2:4,:)')*0.0001; %units em (m/s)^2
    var_gyr = var(Agyr.ans(2:4,:)');
end