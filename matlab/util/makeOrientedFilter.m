function filt = makeOrientedFilter(f,phi, orientation, type, varargin)

    n = varargin{1};
    [x,y,z] = meshgrid(-n:n, -n:n, -n:n);
    r = sqrt(x.^2+y.^2+z.^2);
    
    orientation = orientation / norm(orientation, 2);
    
    anglePrerotation = acos(sum(orientation(:) .* [0;0;1]));
    axis = cross(orientation, [0,0,1]);
    axis = axis/norm(axis);
    rotMat = getDirCosMat(axis,-anglePrerotation);

    coordinates = (horzcat(x(:), y(:), z(:)))';
    rotCoordinates = (rotMat * coordinates)';
    xrot = reshape(rotCoordinates(:,1), size(x));
    yrot = reshape(rotCoordinates(:,2), size(y));
    zrot = reshape(rotCoordinates(:,3), size(z));
    
    [u,v,w] = project2Sphere(xrot, yrot, zrot);
    
    
    angles = acos(w);
    a = reshape(angles, [numel(angles(:)) 1]);
    values = interp1(phi, f, a); % should this use a? test this
    sphericalVals = reshape(values, size(angles));

    g = getRadialFunction3d(type, varargin{:});
    filt = g .* sphericalVals;
end

