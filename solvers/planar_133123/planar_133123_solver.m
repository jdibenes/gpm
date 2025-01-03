
function [R, t, flag] = planar_133123_solver(P11, P12, P21, p22, twp_options)
cof = [sum(p22.^2), -2*dot(p22,P21), sum(P21.^2)-sum((P12-P11).^2)];
sol = roots(cof);
if (imag(sol(1)) < 1e-3), N = 2; else, N = 0; end
sol = real(sol);
R = zeros(3, 3, N);
t = zeros(3, N);
for n = 1:N
    [R(:, :, n), t(:, n)] = twp_solver(P11, P12, P21, sol(n)*p22, twp_options);
end

if (N <= 0)
    flag = -1;
else
    flag = 0;
end
end
