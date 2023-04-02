% compute forward kinematics from twist

function baseToFlange = forwardKinematics(kesi, theta, tool0)

baseToFlange = tool0;
for i = length(theta) : -1 : 1
    baseToFlange = exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), theta(i)) * baseToFlange;
end
end