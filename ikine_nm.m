function [ q ] = ikine_nm(T_d,q0)
%ikine_nm 函数摘要
%  输入T_d为末端相对于基坐标的变换矩阵，矩阵大小6*6；
%  输入q0为初始给定角，单位为弧度，矩阵大小1*6；
%  输出q为关节角的解，单位为弧度，矩阵大小1*6；
%  说明：逆运动求解。利用牛顿迭代法，求解得到一组与初始角相关的解。
%  error为迭代精度，指当前角度的变换矩阵T06c与目标变换矩阵T_d的差值。

error = 1;       %迭代精度，可调

q = q0;  %初始角度

for i = 1:1:1000  %循环迭代1000次
    %计算位置误差
    %首先通过正运动学计算出当前位置
    %每个关节相对于前一关节坐标系的变换矩阵
    T06c = fknie_4dof(q);
        
    f1 = T06c(1, 4) - T_d(1, 4);
    f2 = T06c(2, 4) - T_d(2, 4);
    f3 = T06c(3, 4) - T_d(3, 4);
    f4 = T06c(2, 3) - T_d(2, 3);
    f5 = T06c(3, 3) - T_d(3, 3);
    f6 = T06c(3, 2) - T_d(3, 2);
    f7 = T06c(1, 1) - T_d(1, 1);
    f8 = T06c(1, 2) - T_d(1, 2);
    f9 = T06c(1, 3) - T_d(1, 3);
    f10 = T06c(2, 1) - T_d(2, 1);
    f11 = T06c(2, 2) - T_d(2, 2);
    f12 = T06c(3, 1) - T_d(3, 1);
    F = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12];
%      tic
    df = NewtonMethod_FDiff(T_d,q);
%      toc
    q = q-((df'*df)^(-1)*(df')*F)';
%     q = q-(inv(df'*df+1.4*diag(ones(1,6)))*(df')*F)';
    q = mod(q, 2*pi); 
    
    for j = 1:6
        if abs(q(j)) > pi
            q(j) = q(j) - 2*pi;
        end
    end
    q_error = norm(F);
    if q_error < error
        disp('求解完成')
        break
    end
%     if q_error < 0.00001
%         disp('求解完成')
%         break
%     end
end
    
end

