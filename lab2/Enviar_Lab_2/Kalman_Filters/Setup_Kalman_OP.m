% Setting up of the Kalman Filter described in P. Batista paper
%
%   x(1:3) = omega_meas
%   x(4:6) = acc_meas
%
function Kalman = Setup_Kalman_OP(q, r)
Kalman.C = [eye(3), zeros(3,3)];

Kalman.Q = diag(q);

Kalman.R = diag(r);

