A = 0;
B = 1;
C = 1;
D = 0;
Plant = ss(A,B,C,D);

sys = Plant*[1 1];

Q = 2*10^(-2);
R = 1*10^(-2);
N = 0;

[~,L,~] = kalman(sys, Q, R, N);