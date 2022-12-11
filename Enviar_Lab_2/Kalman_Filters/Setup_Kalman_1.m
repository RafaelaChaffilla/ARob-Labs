% Setting up of the Kalman Filter 1 (without bias)
%
%   x    = theta
%
%   u(1) = omega_meas
%   u(2) = noise(omega_meas)
%
function Kalman = Setup_Kalman_1(q, r, sampleTime)

% Continuous Space State function
A = 0;
B = [1,1];
C = 1;
D = [0,0];

Kalman.sys = ss(A,B,C,D);

Q = diag(q);
R = diag(r);
N = 0;

%Obtain Gain in Continuous System
[~,Kalman.contL,~] = kalman(Kalman.sys, Q, R, N);

% Continuous to Discrete Space Space
Kalman.sys = c2d(Kalman.sys,sampleTime);

% Kalman Space State calculation
[Kalman.sys,Kalman.L,~] = kalman(Kalman.sys, Q, R, N);