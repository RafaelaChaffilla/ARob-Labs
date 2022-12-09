%% Filter Validations
%
%   Duarte Ferro Lopes  , n 95783
%   Jose Medeiro        , n 95811
%   Rafaela Chaffilla   , n 95840
%
%% Initialization
clear all;
close all;
addpath '../Kalman_Filters'
sampleTime  = 0.05;
%% Setup of filters
% 1st filter - no bias

Q           = 2*10^(-2);
R           = 1*10^(-2);
Kalman_1    = Setup_Kalman_1(Q, R, sampleTime);

disp(['1st Kalman Filter gains are L_1 = ' num2str(Kalman_1.L)])
% 2nd filter - bias

Q           = [8*10^(-3),1.5*10^(-3)];
R           = 1*10^(-2);
Kalman_2    = Setup_Kalman_2(Q, R, sampleTime);

disp(['2nd Kalman Filter gains are L_1 = ' num2str(Kalman_2.L(1))...
      '; L_2 = ' num2str(Kalman_2.L(2))]);

% P. Batista filter

Q           = [ ones(1,3)*0.05,...
                ones(1,3)*0.01];
R           = [ 0.05, 0.05, 0.05];
Kalman_OP   = Setup_Kalman_OP(Q, R);

%% Tests
tests = 0;
%
% 1st filter
%   filter choise
filter = 1;
%   angle choise
PHI.gain    = 0;
PHI.bias    = 0;
THETA.gain  = 0;
THETA.bias  = 0;
%   bias choise
bias = [0; 0; 0]*pi/180;
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%   bias choise
bias = [-1; 2; 0]*pi/180;
%   simulation
tests = tests +1;
SIM(2) = sim('Validation_Kalman_Filters');

%
% 2nd filter
%   filter choise
filter = 2;
%   angle choise
PHI.gain    = 0;
PHI.bias    = 0;
THETA.gain  = 0;
THETA.bias  = 0;
%   bias choise
bias = [2; -3; 1]*pi/180;
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%
% 3rd filter
%   filter choise
filter = 3;
%   angle choise
PHI.gain    = 7;
PHI.bias    = 3.5;
THETA.gain  = 20;
THETA.bias  = 0;
%   bias choise
bias = [2; -3; 1];
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%
% 4th filter
%   filter choise
filter = 4;
%   angle choise
PHI.gain    = 7;
PHI.bias    = 3.5;
THETA.gain  = 20;
THETA.bias  = 0;
%   bias choise
bias = [2; -3; 1];
%   simulation
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');
%% Plots

for i=1:tests
    if(i<4)
        figure();
        hold on;
        plot(SIM(i).tout, SIM(i).ang_meas.signals.values(:,2), SIM(i).tout, SIM(i).ang_hat.signals.values(:,2), SIM(i).tout, SIM(i).ang.signals.values(:,2)); 
        legend('\theta medido','\theta estimado','\theta real');
        xlabel("Tempo [s]");
        ylabel("\theta [º]");
        hold off;

        figure();
        hold on;
        plot(SIM(i).tout, SIM(i).ang_meas.signals.values(:,1),SIM(i).tout, SIM(i).ang_hat.signals.values(:,1), SIM(i).tout, SIM(i).ang.signals.values(:,1)); 
        legend('\phi medido','\phi estimado','\phi real');
        xlabel("Tempo [s]");
        ylabel("\phi [º]");
        hold off;
    end
    if(i>3)
        % plot do erro dos angulos
        figure();
        hold on;
        for d =1:SIM(i).ang_hat.signals.dimensions
            plot(SIM(i).ang.time, SIM(i).ang_hat.signals.values(:,d) - SIM(i).ang.signals.values(:,d)); 
        end
        legend('\phi','\theta');
        xlabel("Tempo [s]");
        ylabel('e_{ângulos} [º]');
        hold off;
        % plot do erro dos bias
        figure();
        hold on;
        for d =1:SIM(i).bias_real.signals.dimensions
            plot(SIM(i).bias_hat.time, SIM(i).bias_hat.signals.values(:,d) - SIM(i).bias_real.signals.values(:,d)); 
        end
        legend('$b_x$','$b_y$', '$b_z$','Interpreter','latex');
        xlabel("Tempo [s]");
        ylabel('e_{bias} [º/s]');
        hold off;
    end
end

%% bias
%bias em p
figure();
hold on;
plot(SIM(3).bias_hat.time, 5*ones(size(SIM(3).bias_hat.time, 1)),SIM(3).bias_hat.time, SIM(3).bias_hat.signals.values(:,1)); 
legend('bias real','bias estimado');
xlabel("Tempo [s]");
ylabel("bias [º]");
hold off;

%bias em q
figure();
hold on;
plot(SIM(1,3).bias_hat.time, -3*ones(size(SIM(1,3).bias_hat.time, 1)), SIM(1,3).bias_hat.time, SIM(1,3).bias_hat.signals.values(:,2)); 
legend('bias real','bias estimado');
xlabel("Tempo [s]");
ylabel("bias [º]");

% plot do angulo P. Batista
pallete_1 =["#0072BD";...
            "#D95319";...
            "#009DFF";...
            "#FF611E"];
pallete_2 =["#0072BD";...
            "#D95319";...
            "#638E27";...
            "#009DFF";...
            "#FF611E";...
            "#77AC30"];

figure();
hold on;
for i = 4:5
    for d =1:SIM(i).ang_hat.signals.dimensions
        plot(SIM(i).ang.time, SIM(i).ang_hat.signals.values(:,d) - SIM(i).ang.signals.values(:,d),...
            'Color',pallete_1(d)); 
    end
end
legend('$\phi_c$','$\theta_c$',...
       '$\phi_d$','$\theta_d$',...
       'Interpreter','latex');
xlabel("Tempo [s]");
ylabel('e_{ângulos} [º]');
hold off;
% plot do bias P. Batista
figure();
hold on;
for i = 4:5
    for d =1:SIM(i).bias_real.signals.dimensions
        plot(SIM(i).bias_hat.time, SIM(i).bias_hat.signals.values(:,d) - SIM(i).bias_real.signals.values(:,d),...
            'Color',pallete_2(d)); 
    end
end
legend('$b_x,c$','$b_y,c$', '$b_z,c$',...
       '$b_x,d$','$b_y,d$', '$b_z,d$',...
       'Interpreter','latex');
xlabel("Tempo [s]");
ylabel('e_{bias} [º/s]');
hold off;