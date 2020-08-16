# Trajectory-tracking-
by using velocity Jacoby matrix method and a sliding mode controller, making a trajectory tracking of the robotic arm named Innfos-Gluon-6L3, 


## 文件描述：
1.`Jacob_cross_SDH.m`:求解速度雅各比矩阵，建模采用SDH方法  
2.`Trajectory.m`:利用七段S曲线进行轨迹规划，设计等速趋近率的滑模控制器，以速度雅各比矩阵作为被控对象，实现轨迹跟踪  
3.`fknie_4dof.m`:正运动学求解函数  
4.`ikine_nm.m`:逆运动学求解函数，基于牛顿迭代的方法  
&emsp;`NewtonMethod_FDiff.m`:逆运动学求解函数的子函数，用于建立牛顿迭代方程
