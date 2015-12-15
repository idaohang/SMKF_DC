function [x_pred_SMKF, P_pred_SMKF, S_x_pred_SMKF] = SMKF_pred(x_upd_SMKF, P_upd_SMKF, S_x_upd_SMKF, A, G, q)

x_pred_SMKF = A * x_upd_SMKF;
P_pred_SMKF = A * P_upd_SMKF * A' + q * G * G';
S_x_pred_SMKF = A * S_x_upd_SMKF * A';