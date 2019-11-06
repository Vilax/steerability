function [f_even,f_odd,r]=get_steer2d(n,c)

alpha =[-0.3001 0.7612 -0.6723 0.2388 -0.0275];
k=numel(alpha);
alpha_odd=-alpha(1:2:end);
alpha_even=alpha(2:2:end);

mid=n/2+1;
[x,y]=meshgrid(1:n,1:n);
x=x-mid;
y=y-mid;
[theta,r]=cart2pol(x,y);

r_n=1.5*r/c;
r_n=(r_n-1.5)+1.5;
fr=b_spline(r_n);

ftheta_even=zeros(size(theta));
for i=1:numel(alpha_even),
    exponent=k+1-2*i;
    ftheta_even=ftheta_even+alpha_even(i)*(cos(theta)).^exponent;
end

f_even=ftheta_even.*fr;

% figure(1);
% hold off;
% mesh(f_even)
% hold on;
% figure(2);
% contour(f_even,16);

ftheta_odd=zeros(size(theta));
for i=1:numel(alpha_odd),
    exponent=k+2-2*i;
    ftheta_odd=ftheta_odd+alpha_odd(i)*(cos(theta)).^exponent;
end

f_odd=ftheta_odd.*fr;

% figure(1);
% %hold off;
% mesh(f_odd)




