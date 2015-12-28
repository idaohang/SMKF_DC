function [h1, h2] = plot_axial_len(t, axial_len)

figure
x = (axial_len(:, 1) )';
y = (- axial_len(:, 1) )';
h1 = fill( [t, fliplr(t) ], [x, fliplr(y) ], 'y');
grid on
xlabel('Time')
ylabel('Set-membership uncertainty for state 1')
axis square

figure
x = (axial_len(:, 2) )';
y = (- axial_len(:, 2) )';
h2 = fill( [t, fliplr(t) ], [x, fliplr(y) ], 'y');
grid on
xlabel('Time')
ylabel('Set-membership uncertainty for state 2')
axis square