% compute jacobian end from joint angle
% twist
function jacobianBaseMat = jacobianBase(kesi, theta)
num = length(theta);
jacobianBaseMat = nan(6, num);
T = eye(4, 4);
ad_g(T);
for i = 1 : num
    jacobianBaseMat(:, i) = ad_g(T) * kesi(:, i);
    T = T * exp_se3_r(kesi(1 : 3, i), kesi(4 : 6, i), theta(i));
end
end