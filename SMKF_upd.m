function [x_upd_SMKF, P_upd_SMKF, S_x_upd_SMKF] = SMKF_upd(z_pseudo, x_pred_SMKF, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N)

% p = 0.5;
p = get_p(P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N);

K = ( (1 + inv(p) ) * S_x_pred_SMKF * H' + P_pred_SMKF * H') * inv( (1 + inv(p) ) * H * S_x_pred_SMKF * H' + (1 + p) * S_z_pseudo + H * P_pred_SMKF * H' + r / N);
I = eye(size(K, 1) );
x_upd_SMKF = (I - K * H ) * x_pred_SMKF + K * z_pseudo;
P_upd_SMKF = (I - K * H) * P_pred_SMKF * (I - K * H)' + (r / N) * K * K';
S_x_upd_SMKF = (1 + inv(p) ) * (I - K * H) * S_x_pred_SMKF * (I - K * H)' + (1 + p) * K * S_z_pseudo * K';
