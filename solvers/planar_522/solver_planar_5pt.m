
function [R, t, flag] = solver_planar_5pt(P11, P12, P13, P14, P15, P21, P22, P23, P24, P25, use8, fixeig)
P1 = [P11, P12, P13, P14, P15];
P2 = [P21, P22, P23, P24, P25];

q1 = P1 ./ vecnorm(P1, 2, 1);
q2 = P2 ./ vecnorm(P2, 2, 1);

q = [q1(1,:)'.* q2(1,:)', q1(2,:)'.* q2(1,:)', q1(3,:)'.* q2(1,:)', ...
     q1(1,:)'.* q2(2,:)', q1(2,:)'.* q2(2,:)', q1(3,:)'.* q2(2,:)', ...
     q1(1,:)'.* q2(3,:)', q1(2,:)'.* q2(3,:)', q1(3,:)'.* q2(3,:)'];
 
if (~use8)
q = [q; [1, 0, 0, 0, 1, 0, 0, 0, 1]];
nullSpace = null(q); 
X = nullSpace(:,1);
Y = nullSpace(:,2);
Z = nullSpace(:,3);
else
q(:, 5) = q(:, 5) - q(:, 1);
q(:, 9) = q(:, 9) - q(:, 1);
q = q(:, 2:end);
nullSpace = null(q);
X = nullSpace(:,1);
Y = nullSpace(:,2);
Z = nullSpace(:,3);
X = [-X(4)-X(8); X];
Y = [-Y(4)-Y(8); Y]; 
Z = [-Z(4)-Z(8); Z]; 
end

Q = solver_planar_5pt_Q(X, Y, Z);
xys = math_solve_homogeneous(Q);
xys = math_unit(xys);

x = xys(6);
y = xys(9);
z = xys(10);

E = reshape(x*X + y*Y + z*Z, [3,3]).';

if (fixeig)
[V,D] = eig(E);
d = [abs(D(1,1)),abs(D(2,2)),abs(D(3,3))];
dd = [D(1,1), D(2,2), D(3,3)];
[~, I] = sort(d, 'ascend');
Vs = V(:, I);
ds = dd(:, I);
ds(1) = 0;
ds(2) = ds(2);
ds(3) = -ds(2);
Ds = diag(ds);
E = real(Vs*Ds/Vs);
end

[R, t] = cv_E_to_Rt(E);

for index = 1:4
    t(:, index) = t(:, index).'*(P21 - R(:, :, index)*P11) * t(:, index);
end

flag = 0;
end
