close all;
clear;

%% 0.
% xyz_start=[0,-120.02,533.96];       %�켣���(�ؽڽ�Ϊ0)��ĩ�����꣬��λmm��
xyz_start=[0,-120.02,533.96];
xyz_end=[-30,-45,385];     %�켣�յ��ĩ�����ꣻ
T=10;       %��ɹ켣�滮��ʱ�䣻

%% 1.�켣�滮
L=sqrt((xyz_end(1)-xyz_start(1))^2+(xyz_end(2)-xyz_start(2))^2+(xyz_end(3)-xyz_start(3))^2);
dt=T/7;        %ÿ�ε�ʱ�䳤��
v1=L/(16*dt);   %��һ�μ��ٶȹյ�
J=2*v1/(dt*dt); %�Ӽ��ٶ�
amax=dt*J;      %�����ٶ�
v2=v1+dt*amax;  %�ڶ��μ��ٶȹյ�
vmax=v2+v1;     %�������ٶȹյ�

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

vt=vt1+vt2+vt3+vt4+vt5+vt6+vt7;     %��ʱ���ٶ�

S=zeros(1,length(t));       %��ʱ��λ��
for i=2:length(t)
    S(i)=trapz(t(1:i),vt(1:i));
end 

%��ʱ��xyz��λ��
x_s=xyz_start(1)+(xyz_end(1)-xyz_start(1))/L*S;     
y_s=xyz_start(2)+(xyz_end(2)-xyz_start(2))/L*S;
z_s=xyz_start(3)+(xyz_end(3)-xyz_start(3))/L*S;

%��ʱ��xyz����ٶȷ���
v_x=(xyz_end(1)-xyz_start(1))/L*vt;
v_y=(xyz_end(2)-xyz_start(2))/L*vt;
v_z=(xyz_end(3)-xyz_start(3))/L*vt;

figure(1)
subplot(1,3,1)
plot(t,x_s,'k','linewidth',1);
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('x��λ��','FontName','����','FontSize',12);
title('x��','FontName','����','FontSize',12)
axis square;
grid on;

subplot(1,3,2)
plot(t,y_s,'k','linewidth',1);
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('y��λ��','FontName','����','FontSize',12);
title('y��','FontName','����','FontSize',12)
axis square;
grid on;

subplot(1,3,3)
plot(t,z_s,'k','linewidth',1);
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('z��λ��','FontName','����','FontSize',12);
title('z��','FontName','����','FontSize',12)
axis square;
grid on;

% figure(2)
% plot3(x_s,y_s,z_s,'k','linewidth',2);
% xlabel('x��'),ylabel('y��'),zlabel('z��');
% axis square;
% grid on;

%% 2�켣����.
% %ʹ��ָ�������ɵĻ�ģ����
dth = [0; 0; 0; 0; 0; 0];
% th = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
th = [0; 0; 0; 0; 0; 0];

x=[xyz_start';0;0;0];        %ÿʱ�̵�λ��

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
    %λ�ú��ٶȸ���
    xd=[x_s(i);y_s(i);z_s(i);0;0;0];     %����λ��
    dxd=[v_x(i);v_y(i);v_z(i);0;0;0];
    %�����ǰλ��
    q=th(:, i);
%     x(:, i)=[transl(fknie_4dof(q'));0;0;0];    %��ǰʵ��λ��
    %��⵱ǰ�Ƕ��µ��ſɱȾ���
    Jac = Jacob_cross_SDH(q');
    %���
    e(:, i) = xd - x(:,i);      %���
    s = c*e(:, i);      %��ģ��
    v=dxd + (1/c)*ita*sign(s);    %��е�۵�ĩ��ʵ���ٶ�
    de(:, i) = dxd - v;     %����΢��
    dth(:, i) = inv(Jac+lamda.*diag(ones(1,6)))*v;      %�ؽڽǵ�����
    th(:, i + 1) = th(:, i) + dth(:, i)*0.1;    %��һʱ�̵ĹؽڽǶ�
    x(:, i+1) = x(:, i) + v*0.1;    %��е��ĩ��ʵ��λ��
end


figure(2)
plot(t,e(1:3,:),'linewidth',2);
xlabel('ʱ��/s','FontName','����','FontSize',12);
ylabel('���/mm','FontName','����','FontSize',12);
legend('x���λ�ø������','y���λ�ø������','z���λ�ø������')
title('λ�ø����������','FontName','����','FontSize',12)
axis on;
grid on;

figure(3)
subplot(2, 3, 1);
plot(t, th(1, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��1�ؽڽǶȱ仯����','FontName','����','FontSize',12)

subplot(2, 3, 2);
plot(t, th(2, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��2�ؽڽǶȱ仯����','FontName','����','FontSize',12)

subplot(2, 3, 3);
plot(t, th(3, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��3�ؽڽǶȱ仯����','FontName','����','FontSize',12)

subplot(2, 3, 4);
plot(t, th(4, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��4�ؽڽǶȱ仯����','FontName','����','FontSize',12)

subplot(2, 3, 5);
plot(t, th(5, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��5�ؽڽǶȱ仯����','FontName','����','FontSize',12)

subplot(2, 3, 6);
plot(t, th(6, 1:length(th)-1),'LineWidth',1.5,'color','k');
xlabel('ʱ��','FontName','����','FontSize',12);
ylabel('�Ƕ�ֵ/rad','FontName','����','FontSize',12);
set(gca,'XLim',[min(t) max(t)]);
grid on;
title('��6�ؽڽǶȱ仯����','FontName','����','FontSize',12)