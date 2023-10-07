
clear all

% -------------------------------------------------------------------------

model = 'ackermann';
scale = 10;
N = 10000;

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

R_all = zeros([3, 3, N]);
t_all = zeros([3, N]);
r_all = zeros([4, N]);

for k = 1:N
disp([num2str(k) ' / ' num2str(N)]);
[R, t, r] = pose_random_Rt_planar(scale, model);
R_all(:, :, k) = R;
t_all(:, k) = t;
r_all(:, k) = r;
end

save(['data_' model '_' num2str(scale) '_' num2str(N) '.mat'], 'R_all', 't_all', 'r_all', 'model', 'scale', 'N');
