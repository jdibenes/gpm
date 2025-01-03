
function [R, t, flag, kR, thetaR] = twp_solver(PA1, PB1, PA2, PB2, options)
% [R, t, flag, kR, thetaR] = twp_solver(PA1, PB1, PA2, PB2, options)
% Inputs ------------------------------------------------------------------
% PA1:     3D point A in first image  (3x1)
% PB1:     3D point B in first image  (3x1)
% PA2:     3D point A in second image (3x1)
% PB2:     3D point B in second image (3x1)
% options: Use twp_default_options (or twp_debug_options) to generate a
%          default options structure which can be modified before passing
%          it to this function
% Outputs -----------------------------------------------------------------
% R:      3x3 matrix for rotation
% t:      3x1 vector for translation
% flag:   Negative if solver failed
%         Reasons of failure include bad length correspondence and special
%         geometric cases
% kR:     Axis of rotation
% thetaR: Angle of rotation

R = [];
t = [];

[kR, thetaR, flag] = twp_compute_R(PA1, PB1, PA2, PB2, options);
if (flag < 0), return; end

R = math_axisangle_to_R(kR, thetaR);
t = twp_compute_t(PA1, PB1, PA2, PB2, R, options);
end
