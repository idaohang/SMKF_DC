clear
clc

% m = 20; 
% n = 10; 
% p = 4;
% A = randn(m,n); 
% b = randn(m,1);
% C = randn(p,n); 
% d = randn(p,1); 
% e = rand;
% cvx_begin
%     variable x(n)
%     minimize( norm( A * x - b, 2 ) )
%     subject to
%         C * x == d
%         norm( x, Inf ) <= e
% cvx_end

cvx_begin
    variable x(2)
    K = 1 ./ x(1) + x(2) + 1
    minimize(K^2)
    subject to
        x(1) * x(2) == 1
cvx_end