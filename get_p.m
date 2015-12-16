function p = get_p(P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo)

p = fminbnd(@(p) myfun(p, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo), 1e-2, 1e2);