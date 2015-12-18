function K = get_K(P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N)

I = eye(2);
K0 = ones(2, 1);

fun = @(K) -2 * P_pred_SMKF * H' + 2 * K * H * P_pred_SMKF * H' + 2 * K * r / N -2 * S_x_pred_SMKF * H' + 2 * K * H * S_x_pred_SMKF * H' + 2 * K * S_z_pseudo ...
      + 2 * ( ( -S_x_pred_SMKF * H' + K * H * S_x_pred_SMKF * H') / sqrt(trace( (I - K * H) * S_x_pred_SMKF * (I - K * H)' ) ) + K * S_z_pseudo / sqrt(trace(K * S_z_pseudo * K') ) );
  
K = fsolve(fun, K0);

