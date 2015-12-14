clear
clc
close all

% initial parameters
nsteps = 100;
dT = 1;
A = [1, dT; 0, 1];
G = [ (dT^2) / 2; dT];
H = [1, 0];
q = 10;
r = 10;
N = 10;
N1 = 5;
N2 = 5;
S_Z = 1;

% target
X = zeros(nsteps, size(A, 1) );
x0 = [0; 1];
for i = 1 : nsteps
    if i == 1
        X(i, :) = (A * x0 + G * normrnd(0, sqrt(q) ) )';
    else
        X(i, :) = (A * (X(i - 1, :) )' + G * normrnd(0, sqrt(q) ) )';
    end
end

% measurements
Z = zeros(nsteps, N);
for i = 1 : nsteps
    for j = 1 : N
        Z(i, j) = H * (X(i, :) )' + normrnd(0,  sqrt(r) );
    end
end

% quantization Z2
Z2 = Z(:, (N1 + 1) : (N1 + N2 ) );
lower_bound = min(Z2(:) );
upper_bound = max(Z2(:) );
partition = lower_bound : 2 * sqrt(S_Z) : (upper_bound + 2 * sqrt(S_Z) );
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
        assert(abs(Z2(i, j) - partition(index(i, j) ) ) < sqrt(S_Z) )
    end
end