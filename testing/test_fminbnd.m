clear
clc

P_pred_SMKF = [11.5, 10; 10, 15];
S_x_pred_SMKF = [5, 3; 3, 3];
S_z_pseudo = 0.25;

p = fminbnd(@(p) myfun(p, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo), 1e-2, 1e2);