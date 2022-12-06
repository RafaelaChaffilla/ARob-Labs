% Setting up of the Kalman Filter (with bias)
%
%   x(1) = theta
%   x(2) = bias
%
%   u(1) = omega_meas
%   u(2) = noise(omega_meas)
%   u(3) = noise(bias)
%
function Kalman = Setup_Kalman_2(q, r, sampleTime)

% Continuous Space State function (with noise as input)
A = [0,-1;0,0];
B = [1,1,0;0,0,1];
C = [1,0];
D = [0,0,0];

Kalman.sys = ss(A,B,C,D);
% Continuous to Discrete Space Space
Kalman.sys = c2d(Kalman.sys,sampleTime);

Q = diag(q);
R = diag(r);
N = 0;

% Kalman Space State calculation
[Kalman.sys,Kalman.L,~] = kalman(Kalman.sys, Q, R, N);