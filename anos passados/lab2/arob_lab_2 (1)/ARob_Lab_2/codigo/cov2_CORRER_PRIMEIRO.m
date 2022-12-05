clear;
clc;

%% Get the states and initialize arrays
load('sensor_A.mat')
medicoes_A = sensor;
load('sensor_B.mat')
medicoes_B = sensor;
load('sensor_C.mat')
medicoes_C = sensor;
load('sensor_D.mat')
medicoes_D = sensor;

angulos_g = zeros(size(medicoes_A.time, 1)-1, 3);

load('EXP_A.mat')
estados_A = states;
load('EXP_B.mat')
estados_B = states;
load('EXP_C.mat')
estados_C = states;
load('EXP_D.mat')
estados_D = states;

%% Integration
for i=1:size(medicoes_A.time, 1)-1
    angulos_g(i,:) = trapz(medicoes_A.time(1:i), medicoes_D.signals.values(1:i, 4:6));
end

%% Calculos Intermédios
inicio = 2000;
fim = 5000;

media_A = sum(medicoes_A.signals.values(inicio:fim, :)) / size(medicoes_A.signals.values(inicio:fim, :), 1);
media_B = sum(medicoes_B.signals.values(inicio:fim, :)) / size(medicoes_B.signals.values(inicio:fim, :), 1);
media_C = sum(medicoes_C.signals.values(inicio:fim, :)) / size(medicoes_C.signals.values(inicio:fim, :), 1);
media_D = sum(medicoes_D.signals.values(inicio:fim, :)) / size(medicoes_D.signals.values(inicio:fim, :), 1);

acel_normal = norm(media_A(1:3));
acel_conv = 9.81 / acel_normal;

razao = estados_D.signals.values(inicio:fim, 1) ./ (angulos_g(inicio:fim, 1)); 
n = find(isnan(razao));
razao(n) = 0;
angulos_conv = sum(razao) / size(angulos_g(inicio:fim, 1), 1);


%% Conversoes
medicoes_A_conv(:, 1:3) = acel_conv * medicoes_A.signals.values(:, 1:3); 
medicoes_B_conv(:, 1:3) = acel_conv * medicoes_B.signals.values(:, 1:3); 
medicoes_C_conv(:, 1:3) = acel_conv * medicoes_C.signals.values(:, 1:3); 
medicoes_D_conv(:, 1:3) = acel_conv * medicoes_D.signals.values(:, 1:3);

medicoes_A_conv(:, 4:6) = angulos_conv * medicoes_A.signals.values(:, 4:6);
medicoes_B_conv(:, 4:6) = angulos_conv * medicoes_B.signals.values(:, 4:6);
medicoes_C_conv(:, 4:6) = angulos_conv * medicoes_C.signals.values(:, 4:6);
medicoes_D_conv(:, 4:6) = angulos_conv * medicoes_D.signals.values(:, 4:6);

cov_A_r = cov(medicoes_A_conv(inicio:fim, :));  
media_A_r = sum(medicoes_A_conv(inicio:fim, :)) / size(medicoes_A_conv(inicio:fim, :), 1);
cov_B_r = cov(medicoes_B_conv(inicio:fim, :));  
media_B_r = sum(medicoes_B_conv(inicio:fim, :)) / size(medicoes_B_conv(inicio:fim, :), 1);
cov_C_r = cov(medicoes_C_conv(inicio:fim, :)); 
media_C_r = sum(medicoes_C_conv(inicio:fim, :)) / size(medicoes_C_conv(inicio:fim, :), 1);

%% Angulos do acelerometro
theta_C_est = atan2(medicoes_C.signals.values(inicio:fim, 1), sqrt(medicoes_C.signals.values(inicio:fim, 2).^2 + medicoes_C.signals.values(inicio:fim, 3).^2));
phi_C_est = atan2(-medicoes_C.signals.values(inicio:fim, 2), -medicoes_C.signals.values(inicio:fim, 3));
theta_D_est = atan2(medicoes_D.signals.values(:, 1) - media_C_r(1) / acel_conv, sqrt((medicoes_D.signals.values(:, 2) - media_C_r(2) / acel_conv).^2 + (medicoes_D.signals.values(:, 3)).^2));
phi_D_est = atan2(-medicoes_D.signals.values(:, 2) + media_C_r(2) / acel_conv, -medicoes_D.signals.values(:, 3));

%% Covariancias
cov_pitch_acel = cov(theta_C_est);  
cov_roll_acel = cov(phi_C_est);   
cov_wx = cov(medicoes_C_conv(inicio:fim, 4));
cov_wy = cov(medicoes_C_conv(inicio:fim, 5));
cov_wz = cov(medicoes_C_conv(inicio:fim, 6));
cov_acel = [cov_C_r(1, 1) cov_C_r(2, 2) cov_C_r(3, 3)];

%% Sinais para o Simulink (Ex. 4)
w_ym = [medicoes_D.time medicoes_D_conv(:, 5)];
theta_acel = [medicoes_D.time theta_D_est];
theta_real = [medicoes_D.time estados_D.signals.values(1:6001, 2)] ;


w_xm = [medicoes_D.time medicoes_D_conv(:, 4)];
phi_acel = [medicoes_D.time phi_D_est];
phi_real = [medicoes_D.time estados_D.signals.values(1:6001, 1)];

%% Sinais para o Simulink (Ex. 5)
w_zm = medicoes_D_conv(:, 6);
w = [w_xm(:, 1) w_xm(:, 2) w_ym(:, 2) w_zm];
acel = [theta_acel(:, 1) medicoes_D_conv(:, 1:3) - [media_C_r(1:2) 0]];
acel(1:20, 2:3) = 0;
acel(1:20, 4) = -9.81;
psi_real = [estados_D.signals.values(:, 3); estados_D.signals.values(end, 3)];
euler_real = [phi_real(1:6001, 1) phi_real(1:6001, 2) theta_real(1:6001, 2) psi_real(1:6001)]; 

%% Saves
save('cov_pitch_acel');
save('cov_wy');
save('cov_roll_acel');
save('cov_wx');
save('cov_wz');
save('cov_acel');

save('w_ym')
save('theta_acel');
save('theta_real');

save('w_xm')
save('phi_acel')
save('phi_real')

save('w');
save('acel');
save('euler_real');








