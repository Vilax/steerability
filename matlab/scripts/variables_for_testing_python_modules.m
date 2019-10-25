% this script creates the variables used in testing in python modules
close all; clear all;
steering_3d_filters

save('phi_N4.mat', 'phi')
save('steerAngles_N4.mat', 'steerAngles')
save('bCos_N4.mat', 'bCos')
save('powers_N4.mat', 'powers')
save('orig_filt_N4.mat', 'filt')
save('f_N4.mat', 'f')

% gram matrices
