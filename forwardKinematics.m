% compute forward kinematics from twist

function baseToFlange = forwardKinematics(kesi, theta, tool0)
T = eye(4, 4);
for i = 1 : length(theta)
    T = T * exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), theta(i));
end
baseToFlange = T * tool0;
end