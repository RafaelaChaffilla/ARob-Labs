clear;clc;

%% Sim setup
load('cov_acel');
load('w');
load('euler_real');
load('cov_wx');
load('cov_wy');
load('cov_wz');
Q = diag([100*cov_wx 100*cov_wy 100*cov_wz 0.1*cov_wx 0.1*cov_wy 0.1*cov_wz]);
R = 10*diag([cov_acel(1) cov_acel(2) cov_acel(3)]);
C = [eye(3) zeros(3,3)];
w_0 = [-0.3; -0.3; -0.3];
x_0 = [0;0; -9.81; -0.3; -0.3; -0.3];
A = [0 w_0(3) -w_0(2) 0 x_0(3) -x_0(2);
    -w_0(3) 0 w_0(1) -x_0(3) 0 x_0(1);
    w_0(2) -w_0(1) 0 x_0(2) -x_0(1) 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

var_R = 1;

if var_R == 1
    R = 0.01*R;
end
sim('ricatti', [0 30]);
if var_R == 1
    R = 100*R;
end
load('P');
P = P.Data(:,:,end);

bias = -0.3*ones(3,1);
sim_time = 30;
%% Sim Results
sim('kalman_final', [0 30]);
load('kalman_final_resultados');
load('P_sim');
load('estados_sim')

euler_real_sim = [kalman_final_resultados.Data(1,:); kalman_final_resultados.Data(2,:)];
euler_est_sim = [kalman_final_resultados.Data(3,:); kalman_final_resultados.Data(4,:)];

sse_roll = sum((euler_real_sim(1,:) - euler_est_sim(1,:)).^2);
disp(sse_roll);
sse_pitch = sum((euler_real_sim(2,:) - euler_est_sim(2,:)).^2);
disp(sse_pitch);

figure(1);
plot(kalman_final_resultados.Time', kalman_final_resultados.Data(1,:), kalman_final_resultados.Time',kalman_final_resultados.Data(3,:));
legend('Ângulo de Rolamento real','Ângulo de Rolamento est');
xlabel('Tempo (s)');
ylabel('Ângulo de Rolamento (rad)');

figure(2);
plot(kalman_final_resultados.Time', kalman_final_resultados.Data(2,:), kalman_final_resultados.Time', kalman_final_resultados.Data(4,:));
legend('Ângulo de Picada real','Ângulo de Picada est');
xlabel('Tempo (s)');
ylabel('Pitch (rad)');

figure(3);
plot(kalman_final_resultados.Time', kalman_final_resultados.Data(5,:), kalman_final_resultados.Time', kalman_final_resultados.Data(8,:));
legend('Bias real','Bias est');
xlabel('Tempo (s)');
ylabel('Velocidade angular (rad/s)');

figure(4);
plot(kalman_final_resultados.Time', kalman_final_resultados.Data(6,:), kalman_final_resultados.Time', kalman_final_resultados.Data(9,:));
legend('Bias real','Bias est');
xlabel('Tempo (s)');
ylabel('Velocidade angular (rad/s)');

figure(5);
plot(kalman_final_resultados.Time', kalman_final_resultados.Data(7,:), kalman_final_resultados.Time', kalman_final_resultados.Data(10,:));
legend('Bias real','Bias est');
xlabel('Tempo (s)');
ylabel('Velocidade angular (rad/s)');

figure(6);
estados_sim.Data = estados_sim.Data(:,1:3);
plot(estados_sim)
xlabel('Tempo (s)');
ylabel('Gravidade estimada');
legend('gx','gy','gz')
title('');
print('Gravidade_estimada_inici','-dpng');

figure(7);
plot(P);
xlabel('Tempo (s)');
ylabel('P');
title('');

figure(8)
plot(P_sim);
xlabel('Tempo (s)');
ylabel('P');
title('');