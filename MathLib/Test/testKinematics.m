%TESTKINEMATICS - Unit test for kinematics of iiwa

% main
function tests = testKinematics
clc; close all;
tests = functiontests(localfunctions);
end

% preparation before each test point
function setup(testCase)
testCase.TestData.para = iiwa_config;
end

% test case 1: forward kinematics
function forwardKinematicsTest(testCase)
Config = testCase.TestData.para;

jnt = [-1.2376 1.3495 -2.0080 1.6406 -2.9331 0.2873 1.2626]';
toolExp = [
    0.225241189041861   0.632569577766203   0.741027756594648   0.572225310701854
    0.379618958089047  -0.757443980754261   0.531194938491064  -0.211382122188168
    0.897304771757286   0.161661205314779  -0.410742986888125   0.229479844552782
    0                   0                   0                   1.000000000000000
    ];

toolKesiLearn = forwardKinematics(Config.Tool0, Config.Kesi.Space, jnt);
toolKesiBook = FKinSpace(Config.Tool0, Config.Kesi.Space, jnt);
poseErrorLearn = max(abs(toolKesiLearn - toolExp), [], 'all');
poseErrorBook = max(abs(toolKesiBook - toolExp), [], 'all');

testCase.verifyEqual(poseErrorLearn, 0, 'AbsTol', Config.AbsTol.fair);
testCase.verifyEqual(poseErrorBook, 0, 'AbsTol', Config.AbsTol.fair);
end

% test case 2: jacobianBody and jacobianSpace
function jacobianTest(testCase)
Config = testCase.TestData.para;
numTest = 1e2;
[iiwaJacobianSpace, iiwaJacobianBody, jacobianBodyFromBase, Js, Jb] = deal(nan(6, 7, numTest));
for i = 1 : numTest
    jnt = Config.Range.JntLow + (Config.Range.JntUp - Config.Range.JntLow) .* rand(7, 1);
    toolPos = FKinSpace(Config.Tool0, Config.Kesi.Space, jnt);
    
    iiwaJacobianSpace(:, :, i) = jacobianSpace(Config.Kesi.Space, jnt);
    iiwaJacobianBody(:, :, i) = jacobianBody(Config.Kesi.Body, jnt);
    jacobianBodyFromBase(:, :, i) = Adjoint(inv(toolPos)) * iiwaJacobianSpace(:, :, i);
    Js(:, :, i) = JacobianSpace(Config.Kesi.Space, jnt);
    Jb(:, :, i) = JacobianBody(Config.Kesi.Body, jnt);
end

jacobianSpaceErr = max(abs(iiwaJacobianSpace - Js), [], 'all');
jacobianBodyErr = max(abs(iiwaJacobianBody - Jb), [], 'all');
jacobianBodyErrFromBase = max(abs(jacobianBodyFromBase - Jb), [], 'all');

testCase.verifyEqual(jacobianBodyErr, 0, 'AbsTol', Config.AbsTol.fair);
testCase.verifyEqual(jacobianSpaceErr, 0, 'AbsTol', Config.AbsTol.fair);
testCase.verifyEqual(jacobianBodyErrFromBase, 0, 'AbsTol', Config.AbsTol.fair);
end


% test case 2: jacobianBody and jacobianSpace
function jacobianBaseAnaTest(testCase)
Config = testCase.TestData.para;
numTest = 1e2;
[jacobianCom, jacobian] = deal(nan(6, 7, numTest));
for i = 1 : numTest
    jnt = Config.Range.JntLow + (Config.Range.JntUp - Config.Range.JntLow) .* rand(7, 1);
    jacobianCom(:, :, i) = JacobianAnalytic(Config.Tool0, Config.Kesi.Space, jnt);
    jacobianTmp = iiwa_jacobian_base(jnt);
    jacobianTmp(1 : 3, :) = jacobianTmp(1 : 3, :) / 1e3;
    jacobianTmp = [jacobianTmp(4 : 6, :) ; jacobianTmp(1 : 3, :)];
    jacobian(:, :, i) = jacobianTmp;
end

jacobianErr = max(abs(jacobianCom - jacobian), [], 'all')
testCase.verifyEqual(jacobianErr, 0, 'AbsTol', Config.AbsTol.excellent);
end

% test case 4: jacobianAnalytic
function jacobianAnalyticTest(testCase)
Config = testCase.TestData.para;

sinT = 1;
simT = 1e-5;
numTest = round(sinT / simT);
amp = deg2rad(80) * rand(7, 1);

jnt = nan(7, numTest);
toolPos = nan(4, 4, numTest);
jacobianAnalytic = nan(6, 7, numTest);
t = (1 : numTest) * simT;
for i = 1 : numTest
    jntTmp = amp .* sin(2 * pi / sinT * t(i));
    jnt(:, i) = jntTmp;
    toolPos(:, :, i) = FKinSpace(Config.Tool0, Config.Kesi.Space, jntTmp);
    jacobianAnalytic(:, :, i) = JacobianAnalytic(Config.Tool0, Config.Kesi.Space, jntTmp);
end

%% diff jnt and tool pose
clc;
jntVel = diff(jnt, 1, 2) / simT;
[toolVelInBase, toolVelInBaseCom] = deal(nan(6, numTest - 1));

toolVelInBase(4 : 6, :) = diff(reshape(toolPos(1 : 3, 4, :), 3, []), 1, 2) / simT;
for i = 1 : numTest - 1
    % w = dR * R'
    wTmp = (toolPos(1 : 3, 1 : 3, i + 1) - toolPos(1 : 3, 1 : 3, i)) / simT * toolPos(1 : 3, 1 : 3, i)';
    toolVelInBase(1, i) = wTmp(3, 2);
    toolVelInBase(2, i) = -wTmp(3, 1);
    toolVelInBase(3, i) = wTmp(2, 1);
    
    toolVelInBaseCom(:, i) = jacobianAnalytic(:, :, i) * jntVel(:, i);
end

errorVel = max(abs(toolVelInBaseCom - toolVelInBase), [], 'all');
figure;
for i = 1 : 6
    subplot(2, 3, i);
    hold on; box on;
    plot(toolVelInBaseCom(i, :), '.r');
    plot(toolVelInBase(i, :), '.g');
end
testCase.verifyEqual(errorVel, 0, 'AbsTol', Config.AbsTol.bad);
end