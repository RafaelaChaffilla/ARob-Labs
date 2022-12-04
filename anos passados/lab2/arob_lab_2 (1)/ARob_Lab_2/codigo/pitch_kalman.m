clear;clc;

load('cov_pitch_acel');
load('cov_wy');
load('theta_real');
load('theta_acel');
load('w_ym');

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
Q = cov_wy;
R = 10*cov_pitch_acel;

%Filtro de Kalman [filtro, ganho, covariância]
[kalmf, L, P] = kalman(sys, Q, R);

%% Teste
%Plots
sim('pitch_kalman_sim', [10 30]);
load('pitch_kalman_resultados');

pitch_real_sim = pitch_kalman_resultados.Data(:,1);
pitch_est_sim = pitch_kalman_resultados.Data(:,2);
sse_pitch = sum((pitch_real_sim - pitch_est_sim).^2);
fprintf('sse_pitch = %f\n', sse_pitch);
fprintf('L = %f\n', L);

figure(1);
plot(pitch_kalman_resultados.Time, pitch_kalman_resultados.Data(:, 1), pitch_kalman_resultados.Time, pitch_kalman_resultados.Data(:, 2));
legend('Ângulo de Picada real','Ângulo de Picada estimado');
xlabel('Tempo (s)');
ylabel('Ângulo de Picada (rad)');
