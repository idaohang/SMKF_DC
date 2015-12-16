clear
clc
close all

H = [1, 0];
N = 10;
r = 10;
P_pred_SMKF = [11.5, 10; 10, 15];
S_x_pred_SMKF = [5, 3; 3, 3];
S_z_pseudo = 0.25;

P = 0.1 : 0.1 : 10;
Omega = zeros(length(P), 1);
for i = 1 : length(P)
    p = P(i);
    K = ( (1 + inv(p) ) * S_x_pred_SMKF * H' + P_pred_SMKF * H') * inv( (1 + inv(p) ) * H * S_x_pred_SMKF * H' + (1 + p) * S_z_pseudo + H * P_pred_SMKF * H' + r / N);
    I = eye(size(K, 1) );
    Omega(i) = trace( (I - K * H) * P_pred_SMKF * (I - K * H)' ) + r / N * trace(K * K') + (1 + inv(p) ) * trace( (I - K * H) * S_x_pred_SMKF * (I - K * H)' ) + (1 + p) * trace(K * S_z_pseudo * K');
end

plot(P, Omega)