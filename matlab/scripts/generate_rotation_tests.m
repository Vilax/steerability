% this script makes matrices for testing rotation functions of python

testVol = randn([50, 50, 50]);
save('rotation_test_original.mat', 'testVol', '-v7.3')

axis1 = [1, 3, 4];
axis2 = [-1, 2, 2.1];

rotated1 = rotateFilt3D(testVol, axis1);

rotated2 = rotateFilt3D(testVol, axis2);


save('rotation_test_axis1.mat', 'rotated1', '-v7.3')
save('rotation_test_axi21.mat', 'rotated2', '-v7.3')

% this volume is for testing correct reading in of volumes
% this volume is ones in the front (small in third dimension) and 2 for
% third dimension > 5
easyVol = ones([10, 10, 10]);
easyVol(:,:,6:10) = 2;
save('readin_test_volume.mat', 'easyVol', '-v7.3')

