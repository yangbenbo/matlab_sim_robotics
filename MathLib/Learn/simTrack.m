% simulation of track
%%%%%%%%%%%%%%%%%%%%%%%%%%% change kesi to get error %%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
Config = iiwa_config;

kesiCali = Config.Kesi.Space;
randMat = rand(size(Config.Kesi.Space));

for i = 1 : size(Config.Kesi.Space, 2)
    angleDeg = rand * 1e-2;
    axis = rand(1, 3);
    axis = axis / norm(axis);
    kesiCali(1 : 3, i) = axang2rotm([axis deg2rad(angleDeg)]) * kesiCali(1 : 3, i);
end
kesiCali(4 : 6, :) = kesiCali(4 : 6, :) + rand(3, 7) * 1e-5;


% compare flange error
numTest = 1e4;
for i = 1 : numTest
    jnt = Config.Range.JntLow + (Config.Range.JntUp - Config.Range.JntLow) .* rand(1, 7);
    toolReal = FKinSpace(Config.Tool0, Config.Kesi.Space, jnt);
    toolCali = FKinSpace(Config.Tool0, kesiCali, jnt);
    deltaT = toolReal / toolCali;
    poseErr(i) = norm(deltaT(1 : 3, 4));
    axang = rotm2axang(deltaT(1 : 3, 1 : 3));
    angleErrDeg(i) = rad2deg(axang(4));
end

figure;
subplot(121);
plot(poseErr * 1e3, '.');
subplot(122);
plot(angleErrDeg, '.');

%%
% use kesi calibrate to sim track
clc;
jnt0 = deg2rad([0 30 0 -60 0 90 0]');
baseToToolInit = FKinSpace(Config.Tool0, Config.Kesi.Space, jnt0);
baseToToolDes = baseToToolInit;
baseToToolDes(1 : 3, 4) = baseToToolDes(1 : 3, 4) + rand(3, 1) * 0.01;

angle = deg2rad(rand * 10);
axis = rand(1, 3);
axis = axis / norm(axis);
baseToToolDes(1 : 3, 1 : 3) = baseToToolDes(1 : 3, 1 : 3) * axang2rotm([axis angle]);
deltaT = baseToToolInit \ baseToToolDes;  % transform from current to des: right multiply
deltaPos = deltaT(1 : 3, 4)
axang = rotm2axang(deltaT(1 : 3, 1 : 3));
deltaAngDeg = rad2deg(axang(4))

% use jacobian iterate
maxNum = 1e4;
poseErrThre = 1e-5;
angleErrThre = deg2rad(1e-2);
itI = 0;
jnt = jnt0;
twist = zeros(6, 1);
maxVel = 0.5;
maxAngleVel = deg2rad(1);
simT = 1e-3;
deltaTheta = zeros(7, 1);
[jntArr, poseErrArr, angleErrArr] = deal([]);
while itI < maxNum
    itI = itI + 1;
    jnt = jnt + deltaTheta;
    jacobian = JacobianAnalytic(Config.Tool0, kesiCali, jnt);
    
    BaseToToolCur = FKinSpace(Config.Tool0, Config.Kesi.Space, jnt);
    
    twist(4 : 6) = (baseToToolDes(1 : 3, 4) - BaseToToolCur(1 : 3, 4)) / simT;
    axang = rotm2axang(baseToToolDes(1 : 3, 1 : 3) / BaseToToolCur(1 : 3, 1 : 3));
    twist(1 : 3) = axang(4) * axang(1 : 3) / simT;
    if maxVel < norm(twist(1 : 3))
        twist(1 : 3) = twist(1 : 3) / norm(twist(1 : 3)) * maxVel;
    end
    if maxAngleVel < norm(twist(4 : 6))
        twist(4 : 6) = twist(4 : 6) / norm(twist(4 : 6)) * maxAngleVel;
    end
    deltaTheta = pinv(jacobian) * twist * simT;
    
    poseErr = (norm(baseToToolDes(1 : 3, 4) - BaseToToolCur(1 : 3, 4))) * 1e3;
    angleErr = rad2deg(axang(4));
    
    
    jntArr = [jntArr jnt];
    poseErrArr = [poseErrArr poseErr];
    angleErrArr = [angleErrArr angleErr];
    
    if poseErr < poseErrThre && angleErr < angleErrThre
        poseErr
        angleErr
        disp('finish');
        break;
    end
end

close all;
figure;
subplot(121);
plot(poseErrArr, '-');
subplot(122);
plot(angleErrArr, '-');
figure;
for i = 1 : 7
    subplot(2, 4, i);
    plot(rad2deg(jntArr(i, :)), '-');
end
figure;
for i = 1 : 7
    subplot(2, 4, i);
    jntVel = diff(jntArr(i, :), 1, 2) / simT;
    plot(rad2deg(jntVel), '-');
end