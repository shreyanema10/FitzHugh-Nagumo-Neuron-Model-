% dimensionless model of neuronal excitability V<a
FNmodel(0.3, 0, 0.5, 0.1, 0.1, 0.0, -0.5, 1.5, -0.5, 1.5)

% % dimensionless model of neuronal excitability V>a
% FNmodel(0.6, 0, 0.5, 0.1, 0.1, 0.0, -0.5, 1.5, -.5, 1.5)

% % dimensionless model of neuronal Limit Cycle
% FNmodel(0.6, 0, 0.5, 0.1, 0.1, 0.6, -0.5, 1.5, -0.5, 1.5)

% % dimensionless model of neuronal Depolarization
% FNmodel(0.6, 0, 0.5, 0.1, 0.1, 0.8, -.5, 1.5, -.5, 1.5)

% % dimensionless model of neuronal bistability
% FNmodel(0.6, -0.3, 0.5, 0.8, 0.01, 0.02, -0.5, 1.5, -0.5, 1.5)

function FNmodel(v0, w0, a, r, b, I0, x1, x2, y1, y2) 
% v=Y(1); w=Y(2);
% v0=0.6;w0=0;
% a=0.5;r=0.8;b=0.01;I0=0.02;
% x1 = -.5;
% x2 = 1.5;
% y1 = -.5;
% y2 = 1.5;
Y0=[v0,w0];
t=0:0.1:100;
options=odeset('RelTol',1.e-5);
[T, Y]=ode45(@dydt_FHN,t,Y0,options,a,b,r,I0);
figure(1);clf;
plot(T,Y(:,1),T,Y(:,2)); % time courses of V and w
legend('v(t)','w(t)');
xlabel('\bf{Time}'); ylabel('\bf{v, w}');
title('V(t) vs t and W(t) vs t')
% title(['V(t), W(t) vs t plot for I_{ext}= %d at V(0)= %d', I0, v0]) 

% determine field direction
[v_values, w_values] = meshgrid(x1:0.1:x2, y1:0.1:y2);
vdot= v_values*(a-v_values)*(v_values-1)-w_values+I0;
wdot = b*v_values-r*w_values;

% plot phase plane
figure(2);clf;
hold on;
plot(Y(:,1),Y(:,2)); % V-w phase plane 
%determine and plot the v,w-nullclines
vpts=(x1:.05:x2);
vnullpts=-vpts.*(vpts-a).*(vpts-1)+ I0;
wnullpts= (vpts*b)/r;
plot(vpts,vnullpts,'green',vpts,wnullpts,'red');
hold on;
% plot field direction
h1 = quiver(v_values, w_values, vdot, wdot);
set(h1,'AutoScale','on', 'AutoScaleFactor', 3)
hold on;
% plot starting point
labels = 'Starting point';
plot(v0,w0, '-o')
text(v0,w0,labels, 'VerticalAlignment','top','HorizontalAlignment','left');
% determine fixed point
p = [-1,a+1,-(a+(b/r)),I0];
fp = roots(p);
v = fp(imag(fp) == 0)
w = (v*b)/r
% plot fixed point
hold on;
scatter(v,w,'o')
xlabel('\bf{V}'); ylabel('\bf{W}');
title({'W vs V plot'; 'Phase Plane'})
legend('trajectory', 'v nullcline','w nullcline', 'field direction', 'starting point');
axis([x1 x2 y1 y2]);
end

function dY=dydt_FHN(t,Y,a,b,r,I0)
v=Y(1);
w=Y(2);
dY=zeros(2,1);
dY(1)=-v*(v-a)*(v-1)-w+I0;
dY(2)=b*v-r*w;
end


