function rotatedFilt = rotateFilt3D(filt, targetDir, startDir)
% ROTATEFILT3D  Rotate a 3d filter to new direction
% makes use of imrotate3d; start axis will be north pole if not specified
    
    northpole = [0, 0, 1];
    
    if nargin < 3
        startDir = northpole;
    end
    assert(norm(startDir) > 0);
    startDir = startDir / norm(startDir);
    
    % first compute angle of rotation from 3d in degrees
    targetDir = targetDir / norm(targetDir);
    angle = acos(sum(targetDir(:) .* startDir(:)));
    
    w = cross(startDir, targetDir); % axis of rotation

    rotatedFilt = imrotate3(filt, angle, w, 'linear', 'crop');
    
    
end

