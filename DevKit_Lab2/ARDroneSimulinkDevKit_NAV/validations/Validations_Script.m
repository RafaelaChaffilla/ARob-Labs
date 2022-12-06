clear all;

sampleTime  = 0.05;
%% 1st filter - no bias

A = 0;
B = [1,1];
C = 1;
D = [0,0];

Kalman_1.sys = ss(A,B,C,D);
Kalman_1.sys = c2d(Kalman_1.sys,sampleTime);

Q = 2*10^(-2);
R = 1*10^(-2);
N = 0;

[Kalman_1.sys,Kalman_1.L,~] = kalman(Kalman_1.sys, Q, R, N);

%% 2nd filter - bias

A = [0,-1;0,0];
B = [1,1,0;0,0,1];
C = [1,0];
D = [0,0,0];

Kalman_2.sys = ss(A,B,C,D);
Kalman_2.sys = c2d(Kalman_2.sys,sampleTime);

Q = diag([8*10^(-3),1.5*10^(-3)]);
R = 1*10^(-2);
N = 0;

[Kalman_2.sys,Kalman_2.L,~] = kalman(Kalman_2.sys, Q, R, N);

%% P. Batista filter

Kalman_OP.C = [eye(3), zeros(3,3)];

Kalman_OP.Q = [ eye(3)*0.05 , zeros(3,3)    ;...
                zeros(3,3)  , eye(3)*0.01   ];

Kalman_OP.R = eye(3)*0.05;

SIM.Kalman_OP = sim('Validation_Kalman_OP');