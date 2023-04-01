%% test jacobian
% set joint sine angle to test jacobian
clc; clear; close all;
init_iiwa;
sinT = 1;
simT = 1e-5;
numTest = round(sinT / simT);
amp = deg2rad(80) * rand(7, 1);

jnt = nan(7, numTest);
toolPos = nan(3, numTest);

iiwaJacobianBase = nan(6, 7, numTest);
iiwaJacobianEnd = nan(6, 7, numTest);
t = (1 : numTest) * simT;
for i = 1 : numTest
    jntTmp = amp .* sin(2 * pi / sinT * t(i));
    jnt(:, i) = jntTmp;
    toolPosTmp = forwardKinematics(iiwaTheoryKesi, jntTmp, T_baseToFlangeInit);
    toolPos(:, i) = toolPosTmp(1 : 3, 4);
    iiwaJacobianBase(:, :, i) = jacobianBase(iiwaTheoryKesi, jntTmp);
    iiwaJacobianEnd(:, :, i) = jacobianEnd(iiwaTheoryKesi, jntTmp, T_baseToFlangeInit);
end

jntVel = diff(jnt, 1, 2) / simT;
toolVel = diff(toolPos, 1, 2) / simT;
toolVelCom = nan(size(toolVel));
for i = 1 : numTest - 1
    toolPosTmp = forwardKinematics(iiwaTheoryKesi, jnt(:, i), T_baseToFlangeInit);
    toolVelComTmp = iiwaJacobianEnd(:, :, i) * jntVel(:, i);
    toolVelCom(:, i) = toolPosTmp(1 : 3, 1 : 3) * toolVelComTmp(4 : 6);
end

error = max(abs(toolVelCom - toolVel), [], 'all')

% joint pose
figure(1);
sgtitle('q');
for i = 1 : 7
subplot(2, 4, i);
plot(t, rad2deg(jnt(i, :)), '.');
end

% joint vel
figure(2);
sgtitle('dq');
for i = 1 : 7
subplot(2, 4, i);
plot(t(1 : end - 1), rad2deg(jntVel(i, :)), '.');
end

% end pose
figure(3);
sgtitle('x');
for i = 1 : 3
subplot(2, 4, i);
plot(t, toolPos(i, :), '.');
end

% end vel
figure(4);
subtitle('dx');
for i = 1 : 3
subplot(2, 4, i);
hold on; box on;
plot(t(1 : end - 1), toolVel(i, :), '.g');
plot(t(1 : end - 1), toolVelCom(i, :), '.r');
end



%% compare forward kinematics
clc; clear; close all;
init_iiwa;

numTest = 1e3;
for i = 1 : numTest
    jnt = jntLow + (jntUp - jntLow) .* rand(1, 7);
    toolKesi = forwardKinematics(iiwaTheoryKesi, jnt, T_baseToFlangeInit);
    HomogT = iiwa_forward_kinematics(jnt');
    toolDH = HomogT.T0t;
    toolDH(1 : 3, 4) = toolDH(1 : 3, 4) / 1e3;
    error(i) = max(abs(toolKesi - toolDH), [], 'all');
end

figure;
plot(error, '.');