% P. Batista filter
sampleTime  = 0.05;

Kalman_OP.C = [eye(3), zeros(3,3)];

Kalman_OP.Q = [ eye(3)*0.05 , zeros(3,3)    ;...
                zeros(3,3)  , eye(3)*0.01   ];

Kalman_OP.R = eye(3)*0.05;

SIM.Kalman_OP = sim('Validation_Kalman_OP');