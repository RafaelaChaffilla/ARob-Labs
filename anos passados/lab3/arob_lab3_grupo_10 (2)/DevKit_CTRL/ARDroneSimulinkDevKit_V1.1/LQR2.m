%% LQR
A = [0 1 0;0 0 0; 1 0 0];
B = [0;1;0];
R = 1;
Q = diag([10 1 1]);
C = [1 1 1];
D = 0;
sys = ss(A,B,C,D);
K = lqr(sys,Q,R);

%PSI

A = [0 0;1 0];
B = [1;0];
R = 10;
Q = diag([10 1]);
C = [1 0];
D = 0;
sys = ss(A,B,C,D);
K_psi =  lqr(sys,Q,R);