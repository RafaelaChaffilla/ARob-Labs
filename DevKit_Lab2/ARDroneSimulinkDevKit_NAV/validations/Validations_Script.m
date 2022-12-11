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

Q               = 1*10^(-1);
R               = 5*10^(-1);
Kalman_1    = Setup_Kalman_1(Q, R, sampleTime);

disp(['1st Kalman Filter gains are L_1 = ' num2str(Kalman_1.contL)])
% 2nd filter - bias

Q               = [1*10^(-3),1*10^(-1)];
R               = 5*10^(-1);
Kalman_2    = Setup_Kalman_2(Q, R, sampleTime);

disp(['2nd Kalman Filter gains are L_1 = ' num2str(Kalman_2.contL(1))...
      '; L_2 = ' num2str(Kalman_2.contL(2))]);

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
PHI.gain    = 10;
PHI.bias    = 0;
THETA.gain  = 10;
THETA.bias  = 0;
%   bias choise
bias = [0; 0; 0]; %deg
%   simulation
finalTime = 50;
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%   bias choise
bias = [-1; 2; 0]; %deg
%   simulation
finalTime = 50;
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%
% 2nd filter
%   filter choise
filter = 2;
%   angle choise
PHI.gain    = 10;
PHI.bias    = 0;
THETA.gain  = 10;
THETA.bias  = 0;
%   bias choise
bias = [2; -3; 1]; %deg
%   simulation
finalTime = 50;
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');

%
% 3rd filter
%   filter choise
filter = 3;
%   angle choise
PHI.gain    = 7;
PHI.bias    = 7;
THETA.gain  = 20;
THETA.bias  = 0;
%   bias choise
bias = [2; -3; 1];
%   simulation
finalTime = 100;
tests = tests +1;
SIM(tests) = sim('Validation_Kalman_Filters');
%% Plots
palette.PHI	  =["#0026FF";...
                "#00FFFF"];
palette.THETA =["#FF0000";...
                "#E5BF00"];
palette.BIAS  =["#009DFF";...
                "#FF6A00";...
                "#00CC17"];
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
    if(i == 3) %kalman com bias, print estimação do bias
        % plot do erro dos bias
        figure();
        hold on;
        plot(SIM(i).bias_hat.time, SIM(i).bias_hat.signals.values(:,2),... 
             'Color',palette.BIAS(2)); 
        plot(SIM(i).bias_hat.time, SIM(i).bias_real.signals.values(:,2),... 
             'Color',palette.BIAS(1)); 
        legend('estimado','real');     
        xlabel("Tempo [s]");
        ylabel('bias_q [º/s]');
        hold off;
    end
    if(i>3)
        figure();
        hold on;
        % phi
        plot(SIM(i).ang.time, SIM(i).ang_hat.signals.values(:,1) ...
                             -SIM(i).ang.signals.values(:,1),...
                 'Color',palette.PHI(1));
        % theta
        plot(SIM(i).ang.time, SIM(i).ang_hat.signals.values(:,2) ...
                             -SIM(i).ang.signals.values(:,2),...
                 'Color',palette.THETA(1));
        % detalhes logisticos
        legend('\phi','\theta');
        xlabel("Tempo [s]");
        ylabel('e_{ângulos} [º]');
        hold off;
        % plot do erro dos bias
        figure();
        hold on;
        for d =1:SIM(i).bias_real.signals.dimensions
            plot(SIM(i).bias_hat.time, SIM(i).bias_hat.signals.values(:,d)...
                                      -SIM(i).bias_real.signals.values(:,d),...
                 'Color',palette.BIAS(d)); 
        end
        legend( '$b_x$','$b_y$', '$b_z$',...
                'Interpreter','latex');
        xlabel("Tempo [s]");
        ylabel('e_{bias} [º/s]');
        hold off;
    end
end

% bias
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