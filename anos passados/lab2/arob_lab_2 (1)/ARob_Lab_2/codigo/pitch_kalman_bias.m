clear;clc;

load('cov_pitch_acel');
load('cov_wy');
load('theta_real');
load('theta_acel');
load('w_ym');

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
Q = diag([10000*cov_wy cov_wy]);
R = 10*cov_pitch_acel;

%Filtro de Kalman [filtro, ganho, covariância]
[kalmf, L, P] = kalman(sys, Q, R);
bias = -0.3;

%% Teste
%Plots
sim('pitch_kalman_bias_sim', [10 30]);
load('pitch_kalman_bias_resultados');

pitch_real_sim = pitch_kalman_bias_resultados.Data(:,1);
pitch_est_sim = pitch_kalman_bias_resultados.Data(:,2);
sse_pitch = sum((pitch_real_sim - pitch_est_sim).^2);
disp(sse_pitch);

figure(1);
plot(pitch_kalman_bias_resultados.Time, pitch_kalman_bias_resultados.Data(:, 1), pitch_kalman_bias_resultados.Time, pitch_kalman_bias_resultados.Data(:, 2));
legend('Ângulo de Picada real','Ângulo de Picada estimado');
xlabel('Tempo (s)');
ylabel('Ângulo de Picada (rad)');

figure(2);
plot(pitch_kalman_bias_resultados.Time, pitch_kalman_bias_resultados.Data(:, 3), pitch_kalman_bias_resultados.Time, pitch_kalman_bias_resultados.Data(:, 4));
legend('Bias real','Bias estimado');
xlabel('Tempo (s)');
ylabel('Velocidade angular (rad/s)')