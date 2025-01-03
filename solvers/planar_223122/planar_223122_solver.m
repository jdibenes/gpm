
function [R, t, flag] = planar_223122_solver(P11, P12, p13, p21, p22, p23, swap_px2)
p1 = [P11,P12,p13];
p2 = [p21,p22,p23];
if (~swap_px2)
q_f = solver_solver_prob_pc_relpose_223122(p1, p2);
else
q_f = solver_solver_prob_pc_relpose_223122_b(p1, p2);
end
q_f = cell2mat(q_f);

q_r = real(q_f);
q_i = imag(q_f);

keep = sum(q_i > 1e-3, 1) <= 0;
q_k = q_r(:, keep);
N = size(q_k, 2);
R = zeros(3,3,N);
t = zeros(3,N);

for n = 1:N
    q = q_k(:, n);
    [Rh,d] = math_R_cayley(q(1),q(2),q(3));
    R(:,:,n) = Rh / d;
    v = cross(p2(:,1),p2(:,2));
    At = [cross(R(:,:,n)*p1(:,3),p2(:,3)).';cross(R(:,:,n)*p1(:,2),p2(:,2)).';v.'];
    Bt = [0;0;-v.'*R(:,:,n)*P11];
    t(:,n) = linsolve(At,Bt);
end

if (N <= 0)
    flag = -1;
else
    flag = 0;
end
end
