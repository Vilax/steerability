% build odd polynomial that approximates an even polynomial, a sort of
% parity reverse of 'testQP.m'

n = 100;
dx = 1/n;
x = [0:dx:1]';
k = 4;
v = []; 
for i = k:-1:1
    v = [v x.^i];
end

% H and f matrices
H = v'*v*dx;
f=zeros(k,1);

% No inequality constraints
A = []
b = []

% Equality constraints
Aeq = zeros(2,k);
for i = 1:2:k
    Aeq(1,i) = 1;
end

for i =2:2:k,
    Aeq(2,i)=1;
end
beq = [1,-1];

%Find the minimum
alpha = quadprog(H,f,A,b,Aeq,beq)

%Report the L2 norm
mismatch=sqrt(alpha'*H*alpha)

alpha_even = alpha;
alpha_even(1:2:k)=0;
alpha_odd=alpha;
alpha_odd(2:2:k)=0;
alpha_mono=zeros(size(alpha));
alpha_mono(2)=1;

%Plot the polynomial
y_odd=v*alpha_odd;
y_even=-v*alpha_even;
y_mono=v*alpha_mono;

figure(1);
hold off;
plot(x,y_odd);
hold on;
plot(x,y_even);
hold on
plot(x,y_mono);


figure(10);
plot(x,y_odd-y_even);


%Plot vs angle
theta=acos(x);
figure(2);
hold off;
plot(theta,y_odd);
hold on;
plot(theta,y_even);

%Plot over extended range
x=[-1:dx:1]';
v=[];
for i=k:-1:1,
    v=[v x.^i];
end


y_odd=v*alpha_odd;
y_even=-v*alpha_even;

figure(3);
hold off;
plot(x,y_odd);
hold on;
plot(x,y_even);