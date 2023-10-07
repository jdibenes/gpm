
clear all

%% Settings ---------------------------------------------------------------

% Dataset to use ('general', 'ackermann')
model = 'general';

%% ------------------------------------------------------------------------

% Load poses
scale = 10;
N = 10000;
data = load(['../data/data_' model '_' num2str(scale) '_' num2str(N) '.mat']);

% Select first pose in dataset (out of 10K)
index = 1;

R = data.R_all(:, :, index); % R
t = data.t_all(:, index);    % t

% For comparison with other solvers, all of the generated rotations are
% about the Y axis thus we apply a random change of basis so the rotation 
% axis is not Y (unless we roll the identity)
Rg = pose_random_R();

R = Rg * R * Rg.'; % Ground truth R
t = Rg * t;        % Ground truth t

% 3D points
PA = [-1; 5;  10];
PB = [ 1; 10; 10];

% Generate observations
PA1 = PA;
PB1 = PB;

PA2 = R*PA1 + t;
PB2 = R*PB1 + t;

% Solve
options = twp_debug_options();
[R_estimated, t_estimated] = twp_solver(PA1, PB1, PA2, PB2, options);

% Compute error
disp('Rotation error');
disp(R - R_estimated);

disp('Translation error');
disp(t - t_estimated);
