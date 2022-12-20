%%  =======================================================================
%  Aeronaves Robotizadas - 2022/2023
%  
%  Laboratório 1
% 
%  Grupo 1
%    Hugo Aranha, 95796
%    Mário Cordeiro, 95828
%    Bruno Pedro, 96363
%  ========================================================================
clear all
close all
%% 
load('EXP_A.mat')
sensor_A = sensor_data;
states_A = states;

load('EXP_B.mat')
sensor_B = sensor_data;
states_B = states;

load('EXP_C.mat')
sensor_C = sensor_data;
states_C = states;

load('EXP_D.mat')
sensor_D = sensor_data;
states_D = states;

load('EXP_E.mat')
sensor_E = sensor_data;
states_E = states;

%Adjust data values (sensivity, measurement comes in mV)
% https://www.cdiweb.com/datasheets/te/converting_the_signal_output_of_a_dc_accelerometer_to_acceleration_(g).pdf
 sens_acc = 0.0095; %(EXP A, a_z deve ser 9.81, a partir da diferênça de tensão)
 sens_gir = 0.6697; 

 sensor_A.signals.values(:,4:6) = sens_gir.*sensor_A.signals.values(:,4:6);
 sensor_B.signals.values(:,4:6) = sens_gir.*sensor_B.signals.values(:,4:6);
 sensor_C.signals.values(:,4:6) = sens_gir.*sensor_C.signals.values(:,4:6);
 sensor_D.signals.values(:,4:6) = sens_gir.*sensor_D.signals.values(:,4:6);
 sensor_E.signals.values(:,4:6) = sens_gir.*sensor_E.signals.values(:,4:6);
% 
 sensor_A.signals.values(:,1:3) = sens_acc.*sensor_A.signals.values(:,1:3);
 sensor_B.signals.values(:,1:3) = sens_acc.*sensor_B.signals.values(:,1:3);
 sensor_C.signals.values(:,1:3) = sens_acc.*sensor_C.signals.values(:,1:3);
 sensor_D.signals.values(:,1:3) = sens_acc.*sensor_D.signals.values(:,1:3);
 sensor_E.signals.values(:,1:3) = sens_acc.*sensor_E.signals.values(:,1:3);

%%  Inclinometer_C

roll_raw_C = atan(sensor_C.signals.values(:,2)./sensor_C.signals.values(:,3));
pitch_raw_C = atan(sensor_C.signals.values(:,1)./(sqrt(sensor_C.signals.values(:,2).^2 + sensor_C.signals.values(:,3).^2)));

figure(15)
hold on 
plot(sensor_C.time,roll_raw_C)
plot(states_C.time, states_C.signals.values(:,1))
legend('raw','meas')
xlim([10 30])

close all
%%  Inclinometer_D


roll_raw_D = atan(sensor_D.signals.values(:,2)./sensor_D.signals.values(:,3));
pitch_raw_D = atan(sensor_D.signals.values(:,1)./(sqrt(sensor_D.signals.values(:,2).^2 + sensor_D.signals.values(:,3).^2)));


%%  Inclinometer_E


roll_raw_E = atan(sensor_E.signals.values(:,2)./sensor_E.signals.values(:,3));
pitch_raw_E = atan(sensor_E.signals.values(:,1)./(sqrt(sensor_E.signals.values(:,2).^2 + sensor_E.signals.values(:,3).^2)));

%% KALMAN FILTERS
clc

%dx = Ax+Bu+Gw_1
%y = Cx + Du + Hw_2

A = 0;
B = 1;
C = 1;
D = 0;
G = 1;
H = 0;


Q = 1e-3;
R = 1e2;

sys = ss(A, [B G], C, [D H]);

% Kalman filter
[kalmf, L, P] = kalman(sys, Q, R);

est = estim(sys,L);

state_meas = [sensor_E.time(5:end) sensor_E.signals.values(5:end,5)];
inclinometer = [sensor_E.time(5:end) pitch_raw_E(5:end)];
given_estimation = [states_E.time rad2deg(states_E.signals.values(:,2))];
sim_time = states_E.time(end);

sim('pitch_kalman_sim');


%% Kalman with bias state

%q = wym - by + n1
%dby = 0

A = [0 -1; 0 0];
B = [1;0];
C = [1 0];
D = 0;
H = [0 0];
G = [1 0;0 1];

Q = [1e-4 0;
     0 5e-5];

R = 0.1;


sys = ss(A, [B G], C, [D H]);
[kalmf, L, P] = kalman(sys, Q, R);
est = estim(sys,L);


state_meas = [sensor_E.time(5:end) sensor_E.signals.values(5:end,5)];
inclinometer = [sensor_E.time(5:end) pitch_raw_E(5:end)];
given_estimation = [states_E.time rad2deg(states_E.signals.values(:,2))];
sim_time = states_E.time(end);
sim('pitch_kalman_bias_sim');

%% Implementation of the approach in question 5
T_s = states_E.time(2)-states_E.time(1);% Periodo de amostragem
sim_time = states_E.time(end);
n = sim_time/T_s;


C = [eye(3) zeros(3)];
B = zeros(6); 
D = zeros(3,6);
x = zeros(6,n+1);
x(3,5) = -9.81;

% I'll use the previosly used Q and R matrix, and assume P identity matrix

%cov();

P = (1e-14)*eye(6);

Q_6 = [1e-4 0 0 0 0 0;
       0 1e-4 0 0 0 0;
       0 0 1e-4 0 0 0;
       0 0 0 1e-8 0 0;
       0 0 0 0 1e-8 0;
       0 0 0 0 0 1e-8];

R_3 =[1e-3 0 0;
       0 1e-3 0;
       0 0 1e-3];   
                     % I don't recall the method to derive it from measurements 
                     % (maybe the average of the error in each direction of the gyroscope)

% Implementation of the algorithm as stated in the the slides and in the
% paper mentioned in question 5

L = zeros(6,3);

for p = 5:(n-1)
    % Predição (só se realiza um ciclo,(k = 1))
    w_n = [deg2rad(sensor_E.signals.values(p,4)) deg2rad(sensor_E.signals.values(p,5))...
        deg2rad(sensor_E.signals.values(p,6))];
    
    y_n = [sensor_E.signals.values(p,1) sensor_E.signals.values(p,2)...
        sensor_E.signals.values(p,3)];

    A = [-S(w_n) -S(y_n); zeros(3) zeros(3)];
    sys = ss(A,B,C,D);
    sys_d = c2d(sys,T_s);
    A_d = sys_d.A;

    x(:,p+1) = A_d*x(:,p);

    n_p = 10;
    T_p = T_s/n_p;
    sys_s_p = c2d(sys,T_p);
    A_d_p = sys_s_p.A;
    for ii = 1:n_p
       P = A_d_p*P*(A_d_p') + ((T_p)^2)*Q_6;
    end
    % Correção
    L = P*(C')*(eye(3)\(R_3 + C*P*(C')));
    
    y_n_1 = [sensor_E.signals.values(p,1) sensor_E.signals.values(p,2)...
        sensor_E.signals.values(p,3)];

    x(:,p+1) = x(:,p+1) + L*(y_n_1' - C*x(:,p+1));

    V_subs = (eye(6) - L*C);
    P = V_subs*P*(V_subs') + L*R_3*(L');

end

y_plot = C*x;

roll_y = rad2deg(atan(y_plot(2,:)./y_plot(3,:)));
pitch_y = rad2deg(atan(y_plot(1,:)./(sqrt(y_plot(2,:).^2 + y_plot(3,:).^2))));

figure(12)
hold on
xlim([0,sim_time]);
plot(sensor_E.time(5:end)',(pitch_y(5:end)), sensor_E.time(5:end)', rad2deg(states_E.signals.values(5:end,2)));


hold off

figure(14)
hold on
xlim([0,sim_time]);
plot(sensor_E.time(5:end)',x(5,5:end));


hold off

%%

function M = S(v)
    M = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end