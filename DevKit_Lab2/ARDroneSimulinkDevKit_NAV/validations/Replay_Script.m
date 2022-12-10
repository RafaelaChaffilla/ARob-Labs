%% Replay of Filter
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

%% Loading of files

files = ['../Kalman_Filters/data_for_replay\D.mat';...
         '../Kalman_Filters/data_for_replay\E.mat'];
     
for c=1:size(files,1)
    Temp = load(files(c,:));
    Temp = Temp.ans;
    DATA(c).accs_meas = [Temp(1,:); Temp(2:4,:)]';
    DATA(c).gyro_meas = [Temp(1,:); Temp(5:7,:)]';
    DATA(c).angs_real = [Temp(1,:); Temp(8:9,:)]';
    DATA(c).finalTime = Temp(1,end);
end

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

Q           = [ [1; 1; 1]*0.025;...
                [1; 1; 1]*1*10^(-6)];
R           = [  1; 1; 1]*0.05;
Kalman_OP   = Setup_Kalman_OP(Q, R);

%% Tests
bias        = [0; 0; 1];
seperate    = [1; 1; 0];

% For every replay
for c=1:size(files,1)
    PAR = DATA(c);
    % For every filter
    for d=1:size(bias,1)
        n              	= d + (c-1)*size(bias,1);
        filter          = d;
        SIM.BIAS(n)   	= bias(d);
        SIM.SEPERATE(n) = seperate(d);
        SIM.RESULTS(n)  = sim('Replay_Kalman_Filters');
    end
end

%% Plots
palette.PHI	  =["#0026FF";...
                "#00FFFF"];
palette.THETA =["#FF0000";...
                "#E5BF00"];
palette.BIAS  =["#009DFF";...
                "#FF6A00";...
                "#00CC17"];
palette.RAW    ="#808080";

for i=1:n
    if(SIM.SEPERATE(i) == 1)
        % PHI
        figure();
        hold on;
        % phi Teorico e Estimado
        END = find(SIM.RESULTS(i).tout == 35);
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang_meas.signals.values(1:END,1),...
            'Color',palette.RAW);
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang.signals.values(1:END,1),...
            'Color',palette.PHI(1));
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang_hat.signals.values(1:END,1),...
            'Color',palette.PHI(2));
        
        legend( 'medido','fornecido', 'estimado',...
                'Location','northwest');
        axis([SIM.RESULTS(i).tout(1), SIM.RESULTS(i).tout(END)  ,...
             min([SIM.RESULTS(i).ang_hat.signals.values(1:END,1);...
                  SIM.RESULTS(i).ang.signals.values(1:END,1);
                  SIM.RESULTS(i).ang_meas.signals.values(1:END,1)])*1.2  ,...
             max([SIM.RESULTS(i).ang_hat.signals.values(1:END,1);...
                  SIM.RESULTS(i).ang.signals.values(1:END,1);...
                  SIM.RESULTS(i).ang_meas.signals.values(1:END,1)])*1.2  ]);
        xlabel("Tempo [s]");
        ylabel("\phi [º]");
        hold off;
        % THETA
        figure();
        hold on;
        % theta Teorico e Estimado
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang_meas.signals.values(1:END,2),...
            'Color',palette.RAW);
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang.signals.values(1:END,2),...
            'Color',palette.THETA(1));
        plot(SIM.RESULTS(i).tout(1:END), SIM.RESULTS(i).ang_hat.signals.values(1:END,2),...
            'Color',palette.THETA(2));
        % legenda e outras coisas
        legend( 'medido','fornecido', 'estimado',...
                'Location','northwest');
        axis([SIM.RESULTS(i).tout(1), SIM.RESULTS(i).tout(END)  ,...
             min([SIM.RESULTS(i).ang_hat.signals.values(1:END,2);...
                  SIM.RESULTS(i).ang.signals.values(1:END,2);
                  SIM.RESULTS(i).ang_meas.signals.values(1:END,2)])*1.2  ,...
             max([SIM.RESULTS(i).ang_hat.signals.values(1:END,2);...
                  SIM.RESULTS(i).ang.signals.values(1:END,2);...
                  SIM.RESULTS(i).ang_meas.signals.values(:,2)])*1.2  ]);
        xlabel("Tempo [s]");
        ylabel("\theta [º]");
        hold off;
    else
        figure();
        hold on;
        % phi Teorico e Estimado
        plot(SIM.RESULTS(i).tout, SIM.RESULTS(i).ang.signals.values(:,1),...
            'Color',palette.PHI(1));
        plot(SIM.RESULTS(i).tout, SIM.RESULTS(i).ang_hat.signals.values(:,1),...
            'Color',palette.PHI(2));
        % theta Teorico e Estimado
        plot(SIM.RESULTS(i).tout, SIM.RESULTS(i).ang.signals.values(:,2),...
            'Color',palette.THETA(1));
        plot(SIM.RESULTS(i).tout, SIM.RESULTS(i).ang_hat.signals.values(:,2),...
            'Color',palette.THETA(2));
        % legenda e outras coisas
        legend( '\phi fornecido'    , '\phi estimado'   ,...
                '\theta fornecido'  , '\theta estimado' ,...
                'Location'          , 'northwest'       ,...
                'NumColumns'        , 2);
        axis([SIM.RESULTS(i).tout(1), SIM.RESULTS(i).tout(end)  ,...
             min([SIM.RESULTS(i).ang_hat.signals.values(:,1);...
                  SIM.RESULTS(i).ang.signals.values(:,1);...
                  SIM.RESULTS(i).ang_hat.signals.values(:,2);...
                  SIM.RESULTS(i).ang.signals.values(:,2)])*1.2  ,...
             max([SIM.RESULTS(i).ang_hat.signals.values(:,1);...
                  SIM.RESULTS(i).ang.signals.values(:,1);...
                  SIM.RESULTS(i).ang_hat.signals.values(:,2);...
                  SIM.RESULTS(i).ang.signals.values(:,2)])*1.2  ]);
        xlabel("Tempo [s]");
        ylabel("Ângulo [º]");
        hold off;
    end
    if(SIM.BIAS(i) == 1)
        figure();
        hold on;
        for d =1:SIM.RESULTS(i).bias_hat.signals.dimensions
            plot(SIM.RESULTS(i).bias_hat.time, SIM.RESULTS(i).bias_hat.signals.values(:,d),...
                'Color',palette.BIAS(d)); 
        end
        legend( '$b_x$'         , '$b_y$'       , '$b_z$',...
                'Location'      , 'northwest'   ,...
                'Interpreter'   , 'latex'       );
        axis([SIM.RESULTS(i).tout(1), SIM.RESULTS(i).tout(end)      ,...
             min([SIM.RESULTS(i).bias_hat.signals.values(:,1);...
                  SIM.RESULTS(i).bias_hat.signals.values(:,2);...
                  SIM.RESULTS(i).bias_hat.signals.values(:,3)])*1.2 ,...
             max([SIM.RESULTS(i).bias_hat.signals.values(:,1);...
                  SIM.RESULTS(i).bias_hat.signals.values(:,2);...
                  SIM.RESULTS(i).bias_hat.signals.values(:,3)])*1.2]);
        xlabel("Tempo [s]");
        ylabel('bias [º/s]');
        hold off;
    end
end