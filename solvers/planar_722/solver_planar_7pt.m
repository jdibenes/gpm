
function [R, t, flag] = solver_planar_7pt(P11, P12, P13, P14, P15, P16, P17, P21, P22, P23, P24, P25, P26, P27, use8, fixeig)
P1 = [P11, P12, P13, P14, P15, P16, P17];
P2 = [P21, P22, P23, P24, P25, P26, P27];

q1 = P1(1:3, :) ./ P1(3, :);
q2 = P2(1:3, :) ./ P2(3, :);

q = [q1(1,:)'.* q2(1,:)', q1(2,:)'.* q2(1,:)', q1(3,:)'.* q2(1,:)', ...
     q1(1,:)'.* q2(2,:)', q1(2,:)'.* q2(2,:)', q1(3,:)'.* q2(2,:)', ...
     q1(1,:)'.* q2(3,:)', q1(2,:)'.* q2(3,:)', q1(3,:)'.* q2(3,:)'];
  
if (~use8)
q = [q; [1, 0, 0, 0, 1, 0, 0, 0, 1]];
nullSpace = math_solve_homogeneous(q);%null(q);
else
q(:, 5) = q(:, 5) - q(:, 1);
q(:, 9) = q(:, 9) - q(:, 1);
q = q(:, 2:end);
nullSpace = math_solve_homogeneous(q);%null(q); 
nullSpace = [-nullSpace(4)-nullSpace(8); nullSpace];  
end

E = reshape(nullSpace, [3, 3]).';

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
