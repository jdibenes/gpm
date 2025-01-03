
function [kR, thetaR, flag] = twp_compute_R(PA1, PB1, PA2, PB2, options)
th_iad = options.minimum_length;
th_ied = options.maximum_mismatch;
th_tF  = options.threshold_alignment;
th_dt  = options.minimum_displacement;
th_cos = options.threshold_planarity;
th_tan = options.threshold_axis;
rf_F   = options.refine_initial;
rf_R   = options.refine_final;
th_ra  = options.threshold_refinement;

kR     = [];
thetaR = [];

% validation --------------------------------------------------------------

[V1, L1] = math_unit(PA1 - PB1); % assuming A and B are different
[V2, L2] = math_unit(PA2 - PB2); % V1 and V2 cannot be zero

if (L1 < th_iad), flag = -1; return; end
if (L2 < th_iad), flag = -2; return; end

if (L1 ~= L2)
L12 = sort([L1, L2]);
if ((L12(1) / L12(2)) < th_ied), flag = -3; return; end
end

kR     = V1;
thetaR = 0;

% solve for thetaF --------------------------------------------------------

thetaF = acos(math_clamp(dot(V1, V2), -1, 1)); % [0, pi]

% handle degenerate and special cases
% thetaF =  0: default to pure translation
% thetaF = pi: select any vector orthogonal to V1 (and V2) 
if     (      thetaF  < th_tF), flag = 1; return;
elseif ((pi - thetaF) < th_tF), kF = math_solve_homogeneous([V1.'; V2.']); 
else,                           kF = math_unit(cross(V1, V2)); 
end

% optional: refine RF -----------------------------------------------------

kF     = twp_refine_axis( V1, V2, kF,         rf_F(1), th_ra(1));
thetaF = twp_refine_angle(V1, V2, kF, thetaF, rf_F(2), th_ra(2));

% planarity assumption begins here ----------------------------------------

kR     = kF;
thetaR = thetaF;

[UA, LA] = math_unit(PA1 - PA2); % degenerate case already handled so
[UB, LB] = math_unit(PB1 - PB2); % UA and UB cannot both be zero

% if either is zero the axis of rotation is kF
if ((LA < th_dt) || (LB < th_dt)), flag = 2; return; end

% compute R from axis-angle composition -----------------------------------

tF_2 = thetaF / 2;

cF_2 = cos(tF_2);
sF_2 = sin(tF_2);

ccA = -sF_2*dot(UA, kF); % by the planar motion assumption 
ccB = -sF_2*dot(UB, kF); % both are zero or neither is zero

switch ((abs(ccA) < th_cos) + (abs(ccB) < th_cos))
case 1, flag = 3; return; % one zero and one non-zero (numerical)
case 2, flag = 4; return; % both zero means the axis of rotation is kF
end

kG = cross(kF, V1); % kF and V1 are always orthogonal

% CHOOSE (UA, UB)
% any lineal combination of UA and UB will work except zero
if (LA > LB), cc = ccA; cs = cF_2*dot(UA, V1) + sF_2*dot(UA, kG);
else,         cc = ccB; cs = cF_2*dot(UB, V1) + sF_2*dot(UB, kG);
end

[kR, LR] = math_unit(sF_2*cs*kF + cF_2*cc*V1 + sF_2*cc*kG);
thetaR   = 2*acos(math_clamp((cs*cF_2)/sqrt(cc^2 + cs^2), -1, 1));

if (LR < th_tan), flag = 5; return; end % LR cannot be zero (numerical)

% optional: refine R ------------------------------------------------------

kR     = twp_refine_axis( UA, UB, kR,         rf_R(1), th_ra(1));
thetaR = twp_refine_angle(V1, V2, kR, thetaR, rf_R(2), th_ra(2));

flag = 0;
end

function k = twp_refine_axis(v1, v2, k, steps, th)
if (steps < 1), return; end

for n = 1:steps
JkF = [v1.'; v2.'; 2*k.'];
if (rank(JkF, th) < 3), break; end
k = k - (JkF \ [dot(v1, k); dot(v2, k); sum(k.^2) - 1]);
end

k = math_unit(k);
end

function theta = twp_refine_angle(v1, v2, k, theta, steps, th)
if (steps < 1), return; end

KF = math_v3_to_ssm(k);
F1 = dot(v2, KF*v1);
F2 = dot(v2, (KF^2)*v1);

for n = 1:steps
sinF = sin(theta);
cosF = cos(theta);
df   = (cosF*F2 - sinF*F1);
if (abs(df) < th), break; end
theta = theta - ((cosF*F1 + sinF*F2) / df);
end
end
