%TESTKINEMATICS - Unit test for kinematics of iiwa

% main
function tests = testDynamics
clc; close all;
tests = functiontests(localfunctions);
end

% preparation before each test point
function setup(testCase)
testCase.TestData.para = iiwa_config;
end

% test case 1: inverse dynamics
function inverseDynamicsTest(testCase)
Config = testCase.TestData.para;

thetalist = [0.1; 0.1; 0.1];
dthetalist = [0.1; 0.2; 0.3];
ddthetalist = [2; 1.5; 1];
g = [0; 0; -9.8];
Ftip = [1; 1; 1; 1; 1; 1];
M01 = [[1, 0, 0, 0]; 
    [0, 1, 0, 0]; 
    [0, 0, 1, 0.089159]; 
    [0, 0, 0, 1]];
M12 = [[0, 0, 1, 0.28]; 
    [0, 1, 0, 0.13585]; 
    [-1, 0 ,0, 0]; 
    [0, 0, 0, 1]];
M23 = [[1, 0, 0, 0]; 
    [0, 1, 0, -0.1197]; 
    [0, 0, 1, 0.395]; 
    [0, 0, 0, 1]];
M34 = [[1, 0, 0, 0]; 
    [0, 1, 0, 0]; 
    [0, 0, 1, 0.14225]; 
    [0, 0, 0, 1]];
G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
Glist = cat(3, G1, G2, G3);
Mlist = cat(3, M01, M12, M23, M34); 
Slist = [[1; 0; 1;      0; 1;     0], ...
       [0; 1; 0; -0.089; 0;     0], ...
       [0; 1; 0; -0.089; 0; 0.425]];
taulist = InverseDynamics(thetalist, dthetalist, ddthetalist, g, ...
                        Ftip, Mlist, Glist, Slist);

% Output:
taulistExp = [74.6962 -33.0677 -3.2306]';

testCase.verifyEqual(norm(taulistExp - taulist), 0, 'AbsTol', Config.AbsTol.verypoor);
end