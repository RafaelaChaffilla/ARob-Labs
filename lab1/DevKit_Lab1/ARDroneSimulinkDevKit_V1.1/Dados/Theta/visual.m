%% Visualization of the drone

%% Files

file_exp        = 'State_modded_exp_2_5.mat';
file_command    = 'State_modded_exp_2_5_commands.mat';

exp.raw = load(file_exp);
exp.raw = exp.raw.ans;
%% Data type
%
%   time
%   roll
%   pitch
%   yaw
%   u
%   v
%   w
%   X
%   Y
%   Z

% Position
exp.position.time               = exp.raw(1,:)';
exp.position.signals.values     = [exp.raw(8,:); exp.raw(9,:); exp.raw(10,:)]';
exp.position.signals.dimentions = 3;

% Angle

exp.angle.time                  = exp.raw(1,:)';
exp.angle.signals.values        = [exp.raw(2,:); exp.raw(3,:); exp.raw(4,:)]';
exp.angle.signals.dimentions    = 3;

sim('visual_simulation');