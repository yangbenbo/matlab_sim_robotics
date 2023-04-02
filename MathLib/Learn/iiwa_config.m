%IIWA_CONFIG_PARA - Load parameter for iiwa kinematics module
% 
% Syntax: Config = IIWA_CONFIG_PARA
%
% Outputs:
% Config - Configuration parameter for iiwa kinematics
%     Config.Tool - Coordinate of tool(relative to joint 7)
%     Config.Range - Range of joint angle
%     Config.AbsTol - Absolute computational tolerance
%
% Other m-files required: none
% Subfunctions: getTool, getJointLimits, getCalTolerance
% MAT-files required: none
% 
% See also: JACOBIAN_BODY, JACOBIAN_SPACE
% Author: zilin.tang
% Email: zilin.tang@@whu.edu.cn
% March 2023; Last revision: NAN

function Config = iiwa_config
    
    if nargin >= 1
        error('iiwa_config_para:incorrect input parameter number!');
    end
    
    Config.Tool0 = getTool;
    Config.Range = getJointLimits;
    Config.AbsTol = getCalTolerance;
    Config.Kesi = getKesi;
    
end

function Tool0 = getTool()
    
    % tool coordinate relative to joint 7
    Tool0 = [
                1 0 0 0
                0 1 0 0
                0 0 1 1.306
                0 0 0 1
                ];
end

function Range = getJointLimits
    
    % range of joint angle
    JntUp = deg2rad([170 120 170 120 170 120 175])';
    JntLow = -JntUp;
    Range.JntUp = JntUp;
    Range.JntLow = JntLow;    
end

function AbsTol = getCalTolerance
    
    % absolute computational tolerance
    AbsTol = struct('bad', 1e-3, ...        % resolution of position(mm)
                    'verypoor', 1e-4,...    % resolution of joint angle(deg)
                    'poor', 1e-6, ...       % resolution of joint angle(rad)
                    'fair', 1e-8, ...
                    'good', 1e-10, ...
                    'excellent', 1e-12);
    
end

function Kesi = getKesi
    
    % absolute computational tolerance
    Kesi.Space = [
                0 0 0 0 0 0 0
                0 1 0 -1 0 1 0
                1 0 1 0 1 0 1
                0 -0.36 0 0.78 0 -1.18 0
                0 0 0 0 0 0 0
                0 0 0 0 0 0 0
        ];

    Kesi.Body = Adjoint(inv(getTool)) * Kesi.Space;
end