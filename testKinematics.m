%TESTKINEMATICS - Unit test for kinematics of iiwa

% main
function tests = testKinematics
    tests = functiontests(localfunctions);
    clc;
end

% preparation before each test point
function setup(testCase)
    testCase.TestData.para = iiwa_config;
end

% test case 1: forward kinematics
function forwardKinematicsTest(testCase)
    Config = testCase.TestData.para;
    
    numTest = 1e3;
    for i = 1 : numTest
        jnt = Config.Range.JntLow + (Config.Range.JntUp - Config.Range.JntLow) .* rand(1, 7);
        toolKesi = forwardKinematics(Config.Kesi.Space, jnt, Config.Tool0);
        HomogT = iiwa_forward_kinematics(jnt');
        toolDH = HomogT.T0t;
        toolDH(1 : 3, 4) = toolDH(1 : 3, 4) / 1e3;
        poseError(i) = max(abs(toolKesi - toolDH), [], 'all');
    end

    testCase.verifyEqual(max(poseError), 0, 'AbsTol', Config.AbsTol.fair);
end

% test case 2: jacobianBody and jacobianSpace
function jacobianTest(testCase)
    Config = testCase.TestData.para;
    
    sinT = 1;
    simT = 1e-3;
    numTest = round(sinT / simT);
    amp = deg2rad(80) * rand(7, 1);

    jnt = nan(7, numTest);
    toolPos = nan(3, numTest);

    iiwajacobianSpace = nan(6, 7, numTest);
    iiwajacobianBody = nan(6, 7, numTest);
    jacobianAnalytic = nan(6, 7, numTest);
    Js = nan(6, 7, numTest);
    t = (1 : numTest) * simT;
    transJacobainEndToJacobianAna = zeros(6, 6);
    for i = 1 : numTest
        jntTmp = amp .* sin(2 * pi / sinT * t(i));
        jnt(:, i) = jntTmp;
        toolPosTmp = forwardKinematics(Config.Kesi.Space, jntTmp, Config.Tool0);
        toolPos(:, i) = toolPosTmp(1 : 3, 4);

        iiwajacobianSpace(:, :, i) = jacobianSpace(Config.Kesi.Space, jntTmp);
        iiwajacobianBody(:, :, i) = jacobianBody(Config.Kesi.Body, jntTmp);
        Js(:, :, i) = JacobianSpace(Config.Kesi.Space, jntTmp);
        Jb(:, :, i) = JacobianBody(Config.Kesi.Body, jntTmp);

        transJacobainEndToJacobianAna(1 : 3, 1 : 3) = toolPosTmp(1 : 3, 1 : 3);
        transJacobainEndToJacobianAna(4 : 6, 4 : 6) = toolPosTmp(1 : 3, 1 : 3);
        jacobianAnalytic(:, :, i) = transJacobainEndToJacobianAna * Jb(:, :, i);
    end


    jacobianSpaceErr = max(abs(iiwajacobianSpace - Js), [], 'all');
    jacobianBodyErr = max(abs(iiwajacobianBody - Jb), [], 'all');

    testCase.verifyEqual(jacobianBodyErr, 0, 'AbsTol', Config.AbsTol.fair);
    testCase.verifyEqual(jacobianSpaceErr, 0, 'AbsTol', Config.AbsTol.fair);
end

% test case 3: jacobianAnalytic
function jacobianAnalyticTest(testCase)
    Config = testCase.TestData.para;
    
    sinT = 1;
    simT = 1e-5;
    numTest = round(sinT / simT);
    amp = deg2rad(80) * rand(7, 1);

    jnt = nan(7, numTest);
    toolPos = nan(3, numTest);
    jacobianAnalytic = nan(6, 7, numTest);
    t = (1 : numTest) * simT;
    transJacobainEndToJacobianAna = zeros(6, 6);
    for i = 1 : numTest
        jntTmp = amp .* sin(2 * pi / sinT * t(i));
        jnt(:, i) = jntTmp;
        toolPosTmp = forwardKinematics(Config.Kesi.Space, jntTmp, Config.Tool0);
        toolPos(:, i) = toolPosTmp(1 : 3, 4);

        transJacobainEndToJacobianAna(1 : 3, 1 : 3) = toolPosTmp(1 : 3, 1 : 3);
        transJacobainEndToJacobianAna(4 : 6, 4 : 6) = toolPosTmp(1 : 3, 1 : 3);
        jacobianAnalytic(:, :, i) = transJacobainEndToJacobianAna * JacobianBody(Config.Kesi.Body, jntTmp);
    end

    %% diff jnt and tool pose
    clc;
    jntVel = diff(jnt, 1, 2) / simT;
    
    toolVelInBase = nan(6, numTest - 1);
    toolVelInBaseCom = nan(6, numTest - 1);

    toolVelInBase(4 : 6, :) = diff(toolPos, 1, 2) / simT;
    for i = 1 : numTest - 1
        % w = dR * R'
        toolPosCur = forwardKinematics(Config.Kesi.Space, jnt(:, i), Config.Tool0);
        toolPosNext = forwardKinematics(Config.Kesi.Space, jnt(:, i + 1), Config.Tool0);
        wTmp = (toolPosNext(1 : 3, 1 : 3) - toolPosCur(1 : 3, 1 : 3)) / simT * toolPosCur(1 : 3, 1 : 3)';
        toolVelInBase(1, i) = wTmp(3, 2);
        toolVelInBase(2, i) = -wTmp(3, 1);
        toolVelInBase(3, i) = wTmp(2, 1);

        toolVelInBaseCom(:, i) = jacobianAnalytic(:, :, i) * jntVel(:, i);
    end

    errorVel = max(abs(toolVelInBaseCom - toolVelInBase), [], 'all');
    testCase.verifyEqual(errorVel, 0, 'AbsTol', Config.AbsTol.bad);
end