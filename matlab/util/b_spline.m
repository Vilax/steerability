function f=b_spline(u)
u=min(max(0,u),3);
u1=double(u<1);
u2=double((u>=1)&(u<2));
u3=double(u>=2);

f=0.5*u.^2.*u1;
f=f+0.5*(-3+6*u-2*u.^2).*u2;
f=f+0.5*(3-u).^2.*u3;
