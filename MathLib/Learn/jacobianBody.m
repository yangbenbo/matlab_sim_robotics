function jacobian = jacobianBody(kesi, theta)
% *** VELOCITY KINEMATICS AND STATICS ***
% Takes kesi: [w; v] The joint screw axes in the space frame when the manipulator
%              is at the home position, in the format of a matrix with the
%              screw axes as the columns,
%       theta: A list of joint coordinates. 
% Returns the corresponding space Jacobian (6xn real numbers).
% can use jacobian base to compute
% jacobianEndMat = ad_g(inv(T * tool0)) * jacobainSpace;
% Example Input:
% 
% clear; clc;
% Slist = [[0; 0; 1;   0; 0.2; 0.2], ...
%        [1; 0; 0;   2;   0;   3], ...
%        [0; 1; 0;   0;   2;   1], ...
%        [1; 0; 0; 0.2; 0.3; 0.4]];
% thetalist = [0.2; 1.1; 0.1; 1.2];
% Js = jacobianBody(Slist, thetalist)
% 
% Output:
% Js =
%    -0.0453    0.9950         0    1.0000
%     0.7436    0.0930    0.3624         0
%    -0.6671    0.0362   -0.9320         0
%     2.3259    1.6681    0.5641    0.2000
%    -1.4432    2.9456    1.4331    0.3000
%    -2.0664    1.8288   -1.5887    0.4000

num = length(theta);
T = eye(4, 4);
for i = num : -1 : 1
    jacobian(:, i) = ad_g(T) * kesi(:, i);
    T = T * exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), -theta(i));
end
end