clear all;
close all;

Aacc=load("./dataA/A_acc.mat");
Agyr=load("./dataA/A_gyro.mat");
Bacc=load("./dataB/B_acc.mat");
Bgyr=load("./dataB/B_gyro.mat");
Cacc=load("./dataC/C_acc.mat");
Cgyr=load("./dataC/C_gyro.mat");
PLUSacc=load("./data/PLUS_accelerometer.mat");
PLUSgyr=load("./data/PLUS_gyro.mat");

%% plots
plot_acc_gyr(Aacc, Agyr);
plot_acc_gyr(Bacc.B_acc, Bgyr.B_gyr);
plot_acc_gyr(Cacc.C_acc, Cgyr.C_gyr);
plot_acc_gyr(PLUSacc, PLUSgyr);

%% media e variancia
% remove data from perturbances
Aacc.ans = Aacc.ans(:, 4:end);
Agyr.ans = Agyr.ans(:, 4:end);
Cacc.C_acc.ans = Cacc.C_acc.ans(:, 190:end);
Cgyr.C_gyr.ans = Cgyr.C_gyr.ans(:, 190:end);

[mean_Aacc, var_Aacc, mean_Agyr, var_Agyr] = mean_and_var(Aacc, Agyr);
[mean_Bacc, var_Bacc, mean_Bgyr, var_Bgyr] = mean_and_var(Bacc.B_acc, Bgyr.B_gyr);
[mean_Cacc, var_Cacc, mean_Cgyr, var_Cgyr] = mean_and_var(Cacc.C_acc, Cgyr.C_gyr);
[mean_PLUSacc, var_PLUSacc, mean_PLUSgyr, var_PLUSgyr] = mean_and_var(PLUSacc, PLUSgyr);

%% functions
function plot_acc_gyr(Aacc, Agyr)
    figure();
    hold on;
    for i=2:size(Aacc.ans,1)
        plot(Aacc.ans(1,:), Aacc.ans(i,:).*0.01); %*0.01 para imprimir em m/s
    end
    xlim([Aacc.ans(1,1) Aacc.ans(1,end)]);
    legend('a_x','a_y','a_z');
    xlabel("Tempo [s]");
    ylabel("Acelara��o Linear [m/s]");

    figure();
    hold on;
    for i=2:size(Agyr.ans,1)
        plot(Agyr.ans(1,:), Agyr.ans(i,:));
    end
    xlim([Aacc.ans(1,1) Aacc.ans(1,end)]);
    %legend('\omega_x','\omega_y','\omega_z');
    legend('p','q','r');
    xlabel("Tempo [s]");
    ylabel("Acelara��o Angular [�/s]");

end

function [mean_acc, var_acc, mean_gyr, var_gyr] = mean_and_var(Aacc, Agyr)
    mean_acc = mean(Aacc.ans(2:4,:)')*0.01; %units em m/s
    mean_gyr = mean(Agyr.ans(2:4,:)');
    var_acc = var(Aacc.ans(2:4,:)')*0.0001; %units em (m/s)^2
    var_gyr = var(Agyr.ans(2:4,:)');
end