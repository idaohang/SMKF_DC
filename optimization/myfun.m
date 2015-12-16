function f = myfun(p, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo)

H = [1, 0];
N = 10;
r = 10;

K = ( (1 + inv(p) ) * S_x_pred_SMKF * H' + P_pred_SMKF * H') * inv( (1 + inv(p) ) * H * S_x_pred_SMKF * H' + (1 + p) * S_z_pseudo + H * P_pred_SMKF * H' + r / N);
I = eye(size(K, 1) );
f = trace( (I - K * H) * P_pred_SMKF * (I - K * H)' ) + r / N * trace(K * K') + (1 + inv(p) ) * trace( (I - K * H) * S_x_pred_SMKF * (I - K * H)' ) + (1 + p) * trace(K * S_z_pseudo * K');

