function [h1, h2] = plot_estimation_error(t, x_upd, axial_len, X)

figure
x = (x_upd(:, 1) + axial_len(:, 1) - X(:, 1) )';
y = (x_upd(:, 1) - axial_len(:, 1) - X(:, 1) )';
h1 = fill( [t, fliplr(t) ], [x, fliplr(y) ], 'y');
hold on
plot(t, x_upd(:, 1) - X(:, 1) )
grid on
legend('Set-valued estimator', 'Center of the set-valued estimator')
xlabel('Time')
ylabel('Estimation error for state 1')
axis square

figure
x = (x_upd(:, 2) + axial_len(:, 2) - X(:, 2) )';
y = (x_upd(:, 2) - axial_len(:, 2) - X(:, 2) )';
h2 = fill( [t, fliplr(t) ], [x, fliplr(y) ], 'y');
hold on
plot(t, x_upd(:, 2) - X(:, 2) )
grid on
legend('Set-valued estimator', 'Center of the set-valued estimator')
xlabel('Time')
ylabel('Estimation error for state 2')
axis square

