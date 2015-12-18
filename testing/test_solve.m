clear
clc

H = [1, 0];
N = 10;
r = 10;
I = eye(2);

S_z = 0.25;

syms P S_x  K

eqn = -2 * P * H' + 2 * K * H * P * H' + 2 * K * r / N -2 * S_x * H' + 2 * K * H * S_x * H' + 2 * K * S_z ...
      + 2 * ( ( -S_x * H' + K * H * S_x * H') / sqrt(trace( (I - K * H) * S_x * (I - K * H)' ) ) + K * S_z / sqrt(trace(K * S_z * K') ) ) ...
      == zeros(2, 1);
  
  sol_K = solve(eqn, K)

