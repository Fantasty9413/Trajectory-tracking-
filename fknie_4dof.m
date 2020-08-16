function [ T ] = fknie_4dof( q )
%FKNIE_4DOF 此处显示有关此函数的摘要
%   输入q为弧度,输出T为关节0到关节6的变换矩阵

d=[105.03,0,0,75.66,80.09,44.36];
a=[0,-174.42,-174.42,0,0,0];
alp=[pi/2,0,0,pi/2,-pi/2,0];
offset=[0,-pi/2,0,-pi/2,0,0];
thd=q+offset;

A1=trotz(thd(1))*transl(0,0,d(1))*trotx(alp(1))*transl(a(1),0,0);
A2=trotz(thd(2))*transl(0,0,d(2))*trotx(alp(2))*transl(a(2),0,0);
A3=trotz(thd(3))*transl(0,0,d(3))*trotx(alp(3))*transl(a(3),0,0);
A4=trotz(thd(4))*transl(0,0,d(4))*trotx(alp(4))*transl(a(4),0,0);
A5=trotz(thd(5))*transl(0,0,d(5))*trotx(alp(5))*transl(a(5),0,0);
A6=trotz(thd(6))*transl(0,0,d(6))*trotx(alp(6))*transl(a(6),0,0);

T=A1*A2*A3*A4*A5*A6;

end

