clear;clc;

load('cov_roll_acel');
load('cov_wx');
load('phi_real');
load('phi_acel');
load('w_xm');

%Estados
A = [0 -1;0 0]; %with bias
B = [1;0];
G = [1;1]; %bias has noise
C = [1 0];
D = 0;
H = 0;

%Espaço de estados
sys = ss(A, [B G], C, [D H]);

%matrizes de ajustamento
Q = diag([10000*cov_wx cov_wx]);
R = 10*cov_roll_acel;

%Filtro de Kalman [filtro, ganho, covariância]
[kalmf, L, P] = kalman(sys, Q, R);
bias = -0.3;

%% Teste
%Plots
sim('roll_kalman_bias_sim', [10 30]);
load('roll_kalman_bias_resultados');

roll_real_sim = roll_kalman_bias_resultados.Data(:,1);
roll_est_sim = roll_kalman_bias_resultados.Data(:,2);
sse_roll = sum((roll_real_sim - roll_est_sim).^2);
disp(sse_roll);

figure(1);
plot(roll_kalman_bias_resultados.Time, roll_kalman_bias_resultados.Data(:, 1), roll_kalman_bias_resultados.Time, roll_kalman_bias_resultados.Data(:, 2));
legend('Ângulo de Rolamento real','Ângulo de Rolamento estimado');
xlabel('Tempo (s)');
ylabel('Ângulo de Rolamento (rad)');

figure(2);
plot(roll_kalman_bias_resultados.Time, roll_kalman_bias_resultados.Data(:, 3), roll_kalman_bias_resultados.Time, roll_kalman_bias_resultados.Data(:, 4));
legend('Bias real','Bias estimado');
xlabel('Tempo (s)');
ylabel('Velocidade angular (rad/s)')