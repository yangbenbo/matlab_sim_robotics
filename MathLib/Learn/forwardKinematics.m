% compute forward kinematics from twist

function baseToFlange = forwardKinematics(tool0, kesi, theta)

baseToFlange = tool0;
for i = length(theta) : -1 : 1
    baseToFlange = MatrixExp6(VecTose3(kesi(:, i) * theta(i))) * baseToFlange;
end
end