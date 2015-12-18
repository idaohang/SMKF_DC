function [x_upd_SMKF, P_upd_SMKF, S_x_upd_SMKF, K_upd_SMKF] = SMKF_upd(z_pseudo, x_pred_SMKF, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N)

K_upd_SMKF = get_K(P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N);
p = get_p(S_x_pred_SMKF, S_z_pseudo, K_upd_SMKF, H);

I = eye(size(K_upd_SMKF, 1) );
x_upd_SMKF = (I - K_upd_SMKF * H ) * x_pred_SMKF + K_upd_SMKF * z_pseudo;
P_upd_SMKF = (I - K_upd_SMKF * H) * P_pred_SMKF * (I - K_upd_SMKF * H)' + (r / N) * K_upd_SMKF * K_upd_SMKF';
S_x_upd_SMKF = (1 + inv(p) ) * (I - K_upd_SMKF * H) * S_x_pred_SMKF * (I - K_upd_SMKF * H)' + (1 + p) * K_upd_SMKF * S_z_pseudo * K_upd_SMKF';
