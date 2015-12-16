function p = get_p(P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N)

cvx_begin
    variable p_cvx
%     K = ( (1 + 1 / p_cvx) * S_x_pred_SMKF * H' + P_pred_SMKF * H') * inv( (1 + 1 / p_cvx) * H * S_x_pred_SMKF * H' + (1 + p_cvx) * S_z_pseudo + H * P_pred_SMKF * H' + r / N);
    K = ones(2, 1);
    I = eye(size(K, 1) );
    minimize(  trace( (I - K * H) * P_pred_SMKF * (I - K * H)' ) + (r / N) * trace(K * K') + (1 + 1 / p_cvx) * trace( (I - K * H) * S_x_pred_SMKF * (I - K * H)' ) + (1 + p_cvx) * trace(K * S_z_pseudo * K') )
    subject to
        p_cvx > 0
cvx_end

p = p_cvx;