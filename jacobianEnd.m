% compute jacobian end from joint angle
% twist
function jacobianEndMat = jacobianEnd(kesi, theta, tool0)
num = length(theta);
jacobainBase = nan(6, num);
T = eye(4, 4);
ad_g(T);
for i = 1 : num
    jacobainBase(:, i) = ad_g(T) * kesi(:, i);
    T = T * exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), theta(i));
end

jacobianEndMat = ad_g(inv(T * tool0)) * jacobainBase;
end