%% LQR
%2.2
A = [zeros(3,3) eye(3);zeros(3,6)];
B = [zeros(3,3);eye(3)];
C = [eye(3) zeros(3,3)];
D = zeros(3,3);
sys = ss(A,B,C,D);
R = diag([1 1 1]);
Q = diag([1 1 1 1 1 1]);
K_22 = lqr(sys,Q,R);