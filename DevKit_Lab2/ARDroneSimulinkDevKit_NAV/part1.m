clear all;
close all;

Aacc=load("./data/A_accelerometer.mat");
Agyr=load("./data/A_gyro.mat");
Bacc=load("./data/B_accelerometer.mat");
Bgyr=load("./data/B_gyro.mat");
Cacc=load("./data/C_accelerometer.mat");
Cgyr=load("./data/C_gyro.mat");
PLUSacc=load("./data/PLUS_accelerometer.mat");
PLUSgyr=load("./data/PLUS_gyro.mat");

%% plots
plot_acc_gyr(Aacc, Agyr);
plot_acc_gyr(Bacc, Bgyr);
plot_acc_gyr(Cacc, Cgyr);
plot_acc_gyr(PLUSacc, PLUSgyr);

%% media e variancia

[mean_Aacc, var_Aacc, mean_Agyr, var_Agyr] = mean_and_var(Aacc, Agyr);
[mean_Bacc, var_Bacc, mean_Bgyr, var_Bgyr] = mean_and_var(Bacc, Bgyr);
[mean_Cacc, var_Cacc, mean_Cgyr, var_Cgyr] = mean_and_var(Cacc, Cgyr);
[mean_PLUSacc, var_PLUSacc, mean_PLUSgyr, var_PLUSgyr] = mean_and_var(PLUSacc, PLUSgyr);

%% functions
function plot_acc_gyr(Aacc, Agyr)
    figure();
    hold on;
    for i=2:size(Aacc.ans,1)
        plot(Aacc.ans(1,:), Aacc.ans(i,:));
    end
    legend('a_x','a_y','a_z');

    figure();
    hold on;
    for i=2:size(Agyr.ans,1)
        plot(Agyr.ans(1,:), Agyr.ans(i,:));
    end
    legend('g_x','g_y','g_z');

end

function [mean_acc, var_acc, mean_gyr, var_gyr] = mean_and_var(Aacc, Agyr)
    mean_acc = mean(Aacc.ans(2:4,:)');
    mean_gyr = mean(Agyr.ans(2:4,:)');
    var_acc = var(Aacc.ans(2:4,:)');
    var_gyr = var(Agyr.ans(2:4,:)');
end