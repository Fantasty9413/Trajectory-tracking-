function [ df ] = NewtonMethod_FDiff( T_d,q0 )
%NewtonMethod_FDiff 函数摘要
%   输入T_d为末端相对于基坐标的变换矩阵，矩阵大小6*6;
%   输入q0为逼近角，单位为弧度，矩阵大小1*6;
%   输出df为逼近方程组f的一阶雅各比矩阵，且带入q0求得具体值;
%   说明：用于求解雅各比矩阵df，用于牛顿迭代法求逆运动学。
% tic
d=[105.03,0,0,75.66,80.09,44.36];
a=[0,-174.42,-174.42,0,0,0];
alp=[pi/2,0,0,pi/2,-pi/2,0];
offset=[0,-pi/2,0,-pi/2,0,0];

syms th1 th2 th3 th4 th5 th6;
th=[th1,th2,th3,th4,th5,th6];
thd=offset+th;

T1=trotz(thd(1))*transl(0,0,d(1))*trotx(alp(1))*transl(a(1),0,0);
T2=trotz(thd(2))*transl(0,0,d(2))*trotx(alp(2))*transl(a(2),0,0);
T3=trotz(thd(3))*transl(0,0,d(3))*trotx(alp(3))*transl(a(3),0,0);
T4=trotz(thd(4))*transl(0,0,d(4))*trotx(alp(4))*transl(a(4),0,0);
T5=trotz(thd(5))*transl(0,0,d(5))*trotx(alp(5))*transl(a(5),0,0);
T6=trotz(thd(6))*transl(0,0,d(6))*trotx(alp(6))*transl(a(6),0,0);

T06=T1*T2*T3*T4*T5*T6;
% toc
% % f=T06-T_d;
% % df=zeros(12,6);
% % for i=1:12
% %     for j=1:6
% %        %diff(f(i),th(j));
% %        df(i,j)= diff(f(i),th(j));
% %     end
% % end

f1 = T06(1, 4) - T_d(1, 4);
f2 = T06(2, 4) - T_d(2, 4);
f3 = T06(3, 4) - T_d(3, 4);
f4 = T06(2, 3) - T_d(2, 3);
f5 = T06(3, 3) - T_d(3, 3);
f6 = T06(3, 2) - T_d(3, 2);
f7 = T06(1, 1) - T_d(1, 1);
f8 = T06(1, 2) - T_d(1, 2);
f9 = T06(1, 3) - T_d(1, 3);
f10 = T06(2, 1) - T_d(2, 1);
f11 = T06(2, 2) - T_d(2, 2);
f12 = T06(3, 1) - T_d(3, 1);

% df = [diff(f1, th1), diff(f1, th2), diff(f1, th3), diff(f1, th4), diff(f1, th5), diff(f1, th6);
%     diff(f2, th1), diff(f2, th2), diff(f2, th3), diff(f2, th4), diff(f2, th5), diff(f2, th6);
%     diff(f3, th1), diff(f3, th2), diff(f3, th3), diff(f3, th4), diff(f3, th5), diff(f3, th6);
%     diff(f4, th1), diff(f4, th2), diff(f4, th3), diff(f4, th4), diff(f4, th5), diff(f4, th6);
%     diff(f5, th1), diff(f5, th2), diff(f5, th3), diff(f5, th4), diff(f5, th5), diff(f5, th6);
%     diff(f6, th1), diff(f6, th2), diff(f6, th3), diff(f6, th4), diff(f6, th5), diff(f6, th6);
%     diff(f7, th1), diff(f7, th2), diff(f7, th3), diff(f7, th4), diff(f7, th5), diff(f7, th6);
%     diff(f8, th1), diff(f8, th2), diff(f8, th3), diff(f8, th4), diff(f8, th5), diff(f8, th6);
%     diff(f9, th1), diff(f9, th2), diff(f9, th3), diff(f9, th4), diff(f9, th5), diff(f9, th6);
%     diff(f10, th1), diff(f10, th2), diff(f10, th3), diff(f10, th4), diff(f10, th5), diff(f10, th6);
%     diff(f11, th1), diff(f11, th2), diff(f11, th3), diff(f11, th4), diff(f11, th5), diff(f11, th6);
%     diff(f12, th1), diff(f12, th2), diff(f12, th3), diff(f12, th4), diff(f12, th5), diff(f12, th6)];
% tic
df = jacobian([f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12],th);     %求解F的雅各比矩阵
% toc
df = subs(df,th,q0);
% % df = subs(df, {th1, th2, th3, th4, th5, th6},{q0(1), q0(2), q0(3), q0(4), q0(5), q0(6)});
df = eval(df);

end

