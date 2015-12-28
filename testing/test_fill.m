clear
clc
close all

load log\data_20151228T144704.mat

t = 1 : nsteps;

figure
x = (x_upd(:, 1) + axial_len(:, 1) - X(:, 1) )';
y = (x_upd(:, 1) - axial_len(:, 1) - X(:, 1) )';
fill( [t, fliplr(t) ], [x, fliplr(y) ], 'y')

figure
x1 = (axial_len(:, 1) )';
y1 = (-axial_len(:, 1) )';
fill( [t, fliplr(t) ], [x1, fliplr(y1) ], 'y')

