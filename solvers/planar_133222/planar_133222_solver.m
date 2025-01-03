
function [R, t, flag] = planar_133222_solver(P11, p12, p13, P21, p22, p23)
dP = abs(P11 - P21);
[~, sel] = max(dP);

switch (sel)
    case 1
        cof = @planar_133222_cof_x;
        sol = @planar_133222_sol_x;
        vec = @(P)([0; P(1); P(2)]);
    case 2
        cof = @planar_133222_cof_y;
        sol = @planar_133222_sol_y;
        vec = @(P)([P(1); 0; P(2)]);
    case 3
        cof = @planar_133222_cof_z;
        sol = @planar_133222_sol_z;
        vec = @(P)([P(1); P(2); 0]);
end

Q = cof(P11, p12, p13, P21, p22, p23);

M1 = eq_to_conic(Q(1, :));
M2 = eq_to_conic(Q(2, :));

P = intersectConics(M1, M2);
N = size(P, 2);
if (N > 0), P = P ./ P(3,:); end

R = zeros(3, 3, N);
t = zeros(3, N);

for n = 1:N
    qk = vec(P(:, n));
    qk(sel) = sol(P11, P21, qk);
    [Rh, dh] = math_R_cayley(qk(1), qk(2), qk(3));
    R(:, :, n) = Rh / dh;
    t(:, n) = P21 - R(:, :, n)*P11;
end

if (N <= 0)
    flag = -1;
else
    flag = 0;
end
end

function M = eq_to_conic(eq)
A = eq(1);
B = eq(4);
C = eq(2) / 2;
D = eq(3) / 2;
E = eq(5) / 2;
F = eq(6);

M = [A C D; C B E; D E F]; 
end
