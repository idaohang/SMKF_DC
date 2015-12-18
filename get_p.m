function p = get_p(S_x_pred_SMKF, S_z_pseudo, K_upd_SMKF, H)

I = eye(size(K_upd_SMKF, 1) );

p = sqrt(trace( (I - K_upd_SMKF * H) * S_x_pred_SMKF * (I - K_upd_SMKF * H)' ) ) / sqrt(trace(K_upd_SMKF * S_z_pseudo * K_upd_SMKF') );