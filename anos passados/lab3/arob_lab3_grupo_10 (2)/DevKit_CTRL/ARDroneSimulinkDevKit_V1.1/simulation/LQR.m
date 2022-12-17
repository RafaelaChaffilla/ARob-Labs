A = [zeros(3,3) eye(3);zeros(3,6)];
B = [zeros(3,3);eye(3)];
C = [eye(3) zeros(3,3)];
D = zeros(3,3);
sys = ss(A,B,C,D);

tmp_r = 1;
tmp_q1 = 1;
tmp_q2 = 1;

R = diag([tmp_r tmp_r tmp_r]);
Q = diag([tmp_q1 tmp_q1 tmp_q1 tmp_q2 tmp_q2 tmp_q2]);
K_22 = lqr(sys,Q,R);
sim('ARDroneTTSim22');

traj = zeros(3077, 3);

for i=1:3077
    for j=1:3
        traj(i, j) = Trajectory.signals.values(j, 1, i); 
    end
end

RMSE = sqrt(mean((Position.signals.values(1:3077, :) - traj).^2));
disp(RMSE);

figure;
plot(Position.time(1:3077), Position.signals.values(1:3077, 1), Position.time(1:3077), traj(:, 1));
figure;
plot(Position.time(1:3077), Position.signals.values(1:3077, 2), Position.time(1:3077), traj(:, 2));
figure;
plot(Position.time(1:3077), Position.signals.values(1:3077, 3), Position.time(1:3077), traj(:, 3));