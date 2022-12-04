clear;clc;

load('cov_roll_acel');
load('cov_wx');
load('phi_real');
load('phi_acel');
load('w_xm');

%Estados
A = 0;
B = 1;
C = 1;
D = 0;
G = 1;
H = 0;

%Espaço de estados
sys = ss(A, [B G], C, [D H]);

%matrizes de ajustamento
Q = cov_wx;
R = 10*cov_roll_acel;

%Filtro de Kalman [filtro, ganho, covariância]
[kalmf, L, P] = kalman(sys, Q, R);

%% Teste
%Plots
sim('roll_kalman_sim', [10 30]);
load('roll_kalman_resultados');

roll_real_sim = roll_kalman_resultados.Data(:,1);
roll_est_sim = roll_kalman_resultados.Data(:,2);
sse_roll = sum((roll_real_sim - roll_est_sim).^2);
disp(sse_roll);

figure(1);
plot(roll_kalman_resultados.Time, roll_kalman_resultados.Data(:, 1), roll_kalman_resultados.Time, roll_kalman_resultados.Data(:, 2));
legend('Ângulo de Rolamento real','Ângulo de Rolamento estimado');
xlabel('Tempo (s)');
ylabel('Ângulo de Rolamento (rad)');