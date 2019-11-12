function rotMat = getDirCosMat(axis, angle)
    ux = axis(1);
    uy = axis(2);
    uz = axis(3);

    m11 = cos(angle) + (1 - cos(angle)) * (ux ^ 2);
    m12 = ux*uy*(1 - cos(angle)) - uz * sin(angle);
    m13 = ux*uz*(1 - cos(angle)) + uy * sin(angle);
    m21 = ux*uy*(1 - cos(angle)) + uz * sin(angle);
    m22 = cos(angle) + (1 - cos(angle)) * (uy^2);
    m23 = uy*uz*(1 - cos(angle)) - ux * sin(angle);
    m31 = ux*uz*(1 - cos(angle)) - uy * sin(angle);
    m32 = uy*uz*(1 - cos(angle)) + ux * sin(angle);
    m33 = cos(angle) + (1 - cos(angle)) * (uz^2);

    rotMat = [m11, m12, m13;
                m21, m22, m23;
                m31, m32, m33];

end

