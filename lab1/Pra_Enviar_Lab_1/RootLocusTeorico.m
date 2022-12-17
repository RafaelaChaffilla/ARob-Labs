close all
clear all

% Definição de Parâmetros
k_d = 1;
m = 2;
s = tf('s');

sys = (1 / m)/(s*(s + k_d / m));
figure(1)
rlocus(sys)