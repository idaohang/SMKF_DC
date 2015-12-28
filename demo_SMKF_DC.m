clear
clc
close all

%% initial parameters
nsteps = 50;
dT = 1;
A = [1, dT; 0, 1];
G = [ (dT^2) / 2; dT];
H = [1, 0];
q = 10;
r = 10;
N = 10;
N1 = 5;
N2 = 5;
assert(N == (N1 + N2) )
S_z = 1;

%% target
X = zeros(nsteps, size(A, 1) );
x0 = [0; 1];
for i = 1 : nsteps
    if i == 1
        X(i, :) = (A * x0 + G * normrnd(0, sqrt(q) ) )';
    else
        X(i, :) = (A * (X(i - 1, :) )' + G * normrnd(0, sqrt(q) ) )';
    end
end

%% measurements
Z = zeros(nsteps, N);
for i = 1 : nsteps
    for j = 1 : N
        Z(i, j) = H * (X(i, :) )' + normrnd(0,  sqrt(r) );
    end
end

%% quantization Z2
Z1 = Z(:, 1 : N1);
Z2 = Z(:, (N1 + 1) : (N1 + N2 ) );
lower_bound = min(Z2(:) );
upper_bound = max(Z2(:) );
partition = lower_bound : 2 * sqrt(S_z) : (upper_bound + 2 * sqrt(S_z) );
index = zeros(size(Z2, 1), size(Z2, 2) );
quants = zeros(size(Z2, 1), size(Z2, 2) );
for i = 1 : size(Z2, 1)
    for j = 1 : size(Z2, 2)
        [min_tmp, idx_tmp] = min(abs(Z2(i, j) - partition ) );
        index(i, j) = idx_tmp;
        quants(i, j) = partition(idx_tmp);
    end
end
% testing
for i = 1 : size(Z2, 1)
    for j = 1 : size(Z2, 2)
        assert(abs(Z2(i, j) - quants(i, j) ) < sqrt(S_z) )
    end
end
% calculate Z_pseudo
Z_pseudo = (sum(Z1, 2) + sum(quants, 2) ) / N;
%% set-membership Kalman filter
x_pred = zeros(nsteps, size(X, 2) );
P_pred= zeros(nsteps, size(X, 2), size(X, 2) );
S_x_pred = zeros(nsteps, size(X, 2), size(X, 2) );
x_upd = zeros(nsteps, size(X, 2) );
P_upd = zeros(nsteps, size(X, 2), size(X, 2) );
S_x_upd = zeros(nsteps, size(X, 2), size(X, 2) );
K_upd = zeros(nsteps, 2);
% initial
x0 = [4; 2];
P0 = diag( [4, 5] );
S_x0 = diag( [2, 3] );
% main loop
for i = 1 : nsteps
    % predict
    if i == 1
        x_upd_SMKF = x0;
        P_upd_SMKF = P0;
        S_x_upd_SMKF = S_x0;
    else
        x_upd_SMKF = (x_upd(i - 1, :) )';
        P_upd_SMKF = shiftdim(P_upd(i - 1, :, :) );
        S_x_upd_SMKF = shiftdim(S_x_upd(i - 1, :, :) );
    end 
    [x_pred_SMKF, P_pred_SMKF, S_x_pred_SMKF] = SMKF_pred(x_upd_SMKF, P_upd_SMKF, S_x_upd_SMKF, A, G, q);
    % pseudo-measurement
    z_pseudo = Z_pseudo(i);
    S_z_pseudo = (N2 / N)^2 * S_z;
    % filter
    [x_upd_SMKF, P_upd_SMKF, S_x_upd_SMKF, K_upd_SMKF] = SMKF_upd(z_pseudo, x_pred_SMKF, P_pred_SMKF, S_x_pred_SMKF, S_z_pseudo, H, r, N); 
    % archive
    x_pred(i, :) = x_pred_SMKF';
    P_pred(i, :, :) = P_pred_SMKF;
    S_x_pred(i, :, :) = S_x_pred_SMKF;
    x_upd(i, :) = x_upd_SMKF';
    P_upd(i, :, :) = P_upd_SMKF;
    S_x_upd(i, :, :) = S_x_upd_SMKF;
    K_upd(i, :) = K_upd_SMKF';
end
%% plot
axial_len = zeros(nsteps, 2);
for i = 1 : nsteps
    M = inv(shiftdim(S_x_upd(i, :, :) ) );
    a = M(1, 1);
    b = (M(1, 2) + M(2, 1) ) / 2;
    c = M(2, 2);
    d = - (2 * x_upd(i, 1) * M(1, 1) + x_upd(i, 2) * (M(1, 2) + M(2, 1) ) ) / 2;
    f = - (2 * x_upd(i, 2) * M(2, 2) + x_upd(i, 1) * (M(1, 2) + M(2, 1) ) ) / 2;
    g = (x_upd(i, 1) )^2 * M(1, 1) + (x_upd(i, 2) )^2 * M(2, 2) + x_upd(i, 1) * x_upd(i, 2) * (M(1, 2) + M(2, 1) ) - 1;
    axial_len(i, :) = get_axial_len(a, b, c, d, f, g);
end
t = 1 : nsteps;
[h1, h2] = plot_estimation_error(t, x_upd, axial_len, X);
[h3, h4] = plot_axial_len(t, axial_len);
%% log
logTime = sprintf('%s', datestr(now,30));
if ismac
    dataName = strcat('log/data_', logTime);
else
    dataName = strcat('log\data_', logTime);
end
save(dataName);