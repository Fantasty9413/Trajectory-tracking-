close all;
clear;

%% 0.
% xyz_start=[0,-120.02,533.96];       %轨迹起点(关节角为0)的末端坐标，单位mm；
xyz_start=[0,-120.02,533.96];
xyz_end=[-30,-45,385];     %轨迹终点的末端坐标；
T=10;       %完成轨迹规划的时间；

%% 1.轨迹规划
L=sqrt((xyz_end(1)-xyz_start(1))^2+(xyz_end(2)-xyz_start(2))^2+(xyz_end(3)-xyz_start(3))^2);
dt=T/7;        %每段的时间长度
v1=L/(16*dt);   %第一次加速度拐点
J=2*v1/(dt*dt); %加加速度
amax=dt*J;      %最大加速度
v2=v1+dt*amax;  %第二次加速度拐点
vmax=v2+v1;     %第三次速度拐点

t1 = 1*dt;
t2 = 2*dt;
t3 = 3*dt;
t4 = 4*dt;
t5 = 5*dt;
t6 = 6*dt;
t7 = 7*dt;

t=0:0.1:T;

vt1=1/2*J*t.^2.*(t>=0 & t<t1);
vt2=(v1+amax*(t-t1)).*(t>=t1 & t<t2);
vt3=(vmax-1/2*J*(t3-t).^2).*(t>=t2 & t<t3);
vt4=vmax.*(t>=t3 & t<t4);
vt5=(vmax-1/2*J*(t-t4).^2).*(t>=t4 & t<t5);
vt6=(v2-amax*(t-t5)).*(t>=t5 & t<t6);
vt7=(1/2*J*(t7-t).^2).*(t>=t6 & t<t7);

vt=vt1+vt2+vt3+vt4+vt5+vt6+vt7;     %各时刻速度

S=zeros(1,length(t));       %各时刻位移
for i=2:length(t)
    S(i)=trapz(t(1:i),vt(1:i));
end 

%各时刻xyz的位移
x_s=xyz_start(1)+(xyz_end(1)-xyz_start(1))/L*S;     
y_s=xyz_start(2)+(xyz_end(2)-xyz_start(2))/L*S;
z_s=xyz_start(3)+(xyz_end(3)-xyz_start(3))/L*S;

%各时刻xyz轴的速度分量
v_x=(xyz_end(1)-xyz_start(1))/L*vt;
v_y=(xyz_end(2)-xyz_start(2))/L*vt;
v_z=(xyz_end(3)-xyz_start(3))/L*vt;

figure(1)
subplot(1,3,1)
plot(t,x_s,'k','linewidth',1);
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('x轴位移','FontName','黑体','FontSize',12);
title('x轴','FontName','黑体','FontSize',12)
axis square;
grid on;

subplot(1,3,2)
plot(t,y_s,'k','linewidth',1);
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('y轴位移','FontName','黑体','FontSize',12);
title('y轴','FontName','黑体','FontSize',12)
axis square;
grid on;

subplot(1,3,3)
plot(t,z_s,'k','linewidth',1);
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('z轴位移','FontName','黑体','FontSize',12);
title('z轴','FontName','黑体','FontSize',12)
axis square;
grid on;

% figure(2)
% plot3(x_s,y_s,z_s,'k','linewidth',2);
% xlabel('x轴'),ylabel('y轴'),zlabel('z轴');
% axis square;
% grid on;

%% 2轨迹跟踪.
% %使用指数趋近律的滑模控制
dth = [0; 0; 0; 0; 0; 0];
% th = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
th = [0; 0; 0; 0; 0; 0];

x=[xyz_start';0;0;0];        %每时刻的位姿

% lamda=1;
% lamda=0;
lamda=1;

k = 0.1;
% ita = 0.01;
% ita = 0.0001;
ita = 0.0002;
c = 5;
% c = 0.05;
% alpha = 2;
alpha = 0;

e = [0; 0; 0; 0; 0; 0];
de = zeros(6,1);
for i = 1 : length(t)
    %位置和速度给定
    xd=[x_s(i);y_s(i);z_s(i);0;0;0];     %期望位姿
    dxd=[v_x(i);v_y(i);v_z(i);0;0;0];
    %解出当前位置
    q=th(:, i);
%     x(:, i)=[transl(fknie_4dof(q'));0;0;0];    %当前实际位姿
    %求解当前角度下的雅可比矩阵
    Jac = Jacob_cross_SDH(q');
    %误差
    e(:, i) = xd - x(:,i);      %误差
    s = c*e(:, i);      %滑模面
    v=dxd + (1/c)*ita*sign(s);    %机械臂的末端实际速度
    de(:, i) = dxd - v;     %误差的微分
    dth(:, i) = inv(Jac+lamda.*diag(ones(1,6)))*v;      %关节角的增量
    th(:, i + 1) = th(:, i) + dth(:, i)*0.1;    %下一时刻的关节角度
    x(:, i+1) = x(:, i) + v*0.1;    %机械臂末端实际位姿
end


figure(2)
plot(t,e(1:3,:),'linewidth',2);
xlabel('时间/s','FontName','黑体','FontSize',12);
ylabel('误差/mm','FontName','黑体','FontSize',12);
legend('x轴的位置跟踪误差','y轴的位置跟踪误差','z轴的位置跟踪误差')
title('位置跟踪误差曲线','FontName','黑体','FontSize',12)
axis on;
grid on;

figure(3)
subplot(2, 3, 1);
plot(t, th(1, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第1关节角度变化曲线','FontName','黑体','FontSize',12)

subplot(2, 3, 2);
plot(t, th(2, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第2关节角度变化曲线','FontName','黑体','FontSize',12)

subplot(2, 3, 3);
plot(t, th(3, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第3关节角度变化曲线','FontName','黑体','FontSize',12)

subplot(2, 3, 4);
plot(t, th(4, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第4关节角度变化曲线','FontName','黑体','FontSize',12)

subplot(2, 3, 5);
plot(t, th(5, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第5关节角度变化曲线','FontName','黑体','FontSize',12)

subplot(2, 3, 6);
plot(t, th(6, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('时间','FontName','黑体','FontSize',12);
ylabel('角度值/rad','FontName','黑体','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('第6关节角度变化曲线','FontName','黑体','FontSize',12)