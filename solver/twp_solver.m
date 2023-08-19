
function [R, t, flag, kR, thetaR] = twp_solver(PA1, PB1, PA2, PB2, options)
R = [];
t = [];

[kR, thetaR, flag] = twp_compute_R(PA1, PB1, PA2, PB2, options);
if (flag < 0), return; end

R = math_axisangle_to_R(kR, thetaR);
t = twp_compute_t(PA1, PB1, PA2, PB2, R, options);
end
