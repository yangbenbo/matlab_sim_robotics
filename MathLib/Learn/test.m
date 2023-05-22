% test [r][w][r]w = [w][r][r]w
% because of [w][r]-[r][w] = [[w]r]
clc; clear; close all;
syms a b c x y z r11 r12 r13 r21 r22 r23 r31 r32 r33 real
syms I11 I12 I13 I21 I22 I23 I31 I32 I33 real

r = [a b c]'
w = [x y z]'
R = [r11 r12 r13 
    r21 r22 r23 
    r31 r32 r33
    ]
R = eye(3, 3)
I = [I11 I12 I13 
    I12 I22 I23 
    I13 I23 I33
    ]
r_screw = VecToso3(r);
w_screw = VecToso3(w);

err = (r_screw * w_screw * r_screw - w_screw * r_screw * r_screw) * w;
simplify(err)
% pretty(simplify(err))

r_screw * w_screw - w_screw * r_screw

err = (r_screw * w_screw - w_screw * r_screw) * w;
simplify(err)
% pretty(simplify(err))


%% 验证欧拉运动方程
err2 = (R'* (I * w_screw - w_screw * I) * R - w_screw * R' * I * R) * w
simplify(err2)
% pretty(simplify(err2))