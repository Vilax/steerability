function [u, v, w] = project2Sphere(x,y,z)
    SMALL_CONSTANT = 1e-6;
    r = sqrt(x.^2 + y.^2 + z.^2 + SMALL_CONSTANT);
    u = x ./ r;
    v = y ./ r;
    w = z ./ r;
%     u(r==0)=0;
%     v(r==0)=0;
%     w(r==0)=0;
end

