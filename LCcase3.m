function LCcase3
Y1=[0.8452,0.5000];
Y2=[0.7,0.5500];
Y3=[0.654,1];
Y4=[1.2, 0.8542];
Y5=[0.6,0.5000];
Y6=[1.3,1.2];
a=0.5;r=0.1;b=0.1;I0=0.8;
x1 = -1.5;
x2 = 1.5;
y1 = -1.5;
y2 = 1.5;
t=0:0.1:100;
options=odeset('RelTol',1.e-5);
[T1, Yt1]=ode45(@dydt_FHN,t,Y1,options,a,b,r,I0);
options=odeset('RelTol',1.e-5);
[T2, Yt2]=ode45(@dydt_FHN,t,Y2,options,a,b,r,I0);
options=odeset('RelTol',1.e-5);
[T3, Yt3]=ode45(@dydt_FHN,t,Y3,options,a,b,r,I0);
options=odeset('RelTol',1.e-5);
[T4, Yt4]=ode45(@dydt_FHN,t,Y4,options,a,b,r,I0);
options=odeset('RelTol',1.e-5);
[T5, Yt5]=ode45(@dydt_FHN,t,Y5,options,a,b,r,I0);
options=odeset('RelTol',1.e-5);
[T6, Yt6]=ode45(@dydt_FHN,t,Y6,options,a,b,r,I0);


[v_values, w_values] = meshgrid(x1:0.1:x2, y1:0.1:y2);
vdot= v_values*(a-v_values)*(v_values-1)-w_values+I0;
wdot = b*v_values-r*w_values;
%determine and plot the v,w-nullclines
vpts=(x1:.05:x2);
vnullpts=-vpts.*(vpts-a).*(vpts-1)+ I0;
wnullpts= (vpts*b)/r;
figure(1);clf;
hold on;
plot(vpts,vnullpts,'green',vpts,wnullpts,'red');
hold on;
h1 = quiver(v_values, w_values, vdot, wdot);
set(h1,'AutoScale','on', 'AutoScaleFactor', 3)
hold on;
plot(Yt1(:,1),Yt1(:,2)); % V-w phase plane
hold on;
plot(Yt2(:,1),Yt2(:,2)); % V-w phase plane
hold on;
plot(Yt3(:,1),Yt3(:,2)); % V-w phase plane
hold on;
plot(Yt4(:,1),Yt4(:,2)); % V-w phase plane
hold on;
plot(Yt5(:,1),Yt5(:,2)); % V-w phase plane
hold on;
plot(Yt6(:,1),Yt6(:,2)); % V-w phase plane
hold on;
plot(Y1(1),Y1(2), '-o')
hold on;
plot(Y2(1),Y2(2), '-o')
hold on;
plot(Y3(1),Y3(2), '-o')
hold on;
plot(Y4(1),Y4(2), '-o')
hold on;
plot(Y5(1),Y5(2), '-o')
hold on;
plot(Y6(1),Y6(2), '-o')
% labels = 'Starting point';
% plot(v0,w0, '-o')
% text(v0,w0,labels, 'VerticalAlignment','top','HorizontalAlignment','left');
xlabel('\bf{V}'); ylabel('\bf{W}');
title({'W vs V plot for V(0)>a at I_{ext}=0.8'; 'Phase Plane with Limit Cycle'})
legend('v nullcline','w nullcline', 'field direction', 'limit cycle trajectory');
axis([-.5 x2 0.3 y2]);
end
function dY=dydt_FHN(t,Y,a,b,r,I0)
v=Y(1);
w=Y(2);
dY=zeros(2,1);
dY(1)=-v*(v-a)*(v-1)-w+I0;
dY(2)=b*v-r*w;
end

