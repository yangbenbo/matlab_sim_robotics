function jacobian = jacobianSpace(kesi, theta)
% *** VELOCITY KINEMATICS AND STATICS ***
% Takes kesi: [w; v], The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       kesi: A list of joint coordinates. 
% Returns the corresponding space Jacobian (6xn real numbers).
% Example Input:
% 
% clear; clc;
% kesi = [[0; 0; 1;   0; 0.2; 0.2], ...
%        [1; 0; 0;   2;   0;   3], ...
%        [0; 1; 0;   0;   2;   1], ...
%        [1; 0; 0; 0.2; 0.3; 0.4]];
% theta = [0.2; 1.1; 0.1; 1.2];
% Js = jacobianSpace(kesi, theta)
% 
% Output:
% Js =
%         0    0.9801   -0.0901    0.9575
%         0    0.1987    0.4446    0.2849
%    1.0000         0    0.8912   -0.0453
%         0    1.9522   -2.2164   -0.5116
%    0.2000    0.4365   -2.4371    2.7754
%    0.2000    2.9603    3.2357    2.2251


num = length(theta);
T = eye(4, 4);
for i = 1 : num
    jacobian(:, i) = Adjoint(T) * kesi(:, i);
%     T = T * exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), theta(i));
    T = T * MatrixExp6(VecTose3(kesi(:, i) * theta(i)));
end
end