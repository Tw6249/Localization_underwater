% 测试三维欧式空间中的矩阵的秩
clear;
clc;

% 随机生成3维空间中的点
num_points = 100;
points = rand(num_points,3)*10; % 在20*20*20的空间中随机生成点

% 绘图
% scatter3(points(:,1),points(:,2),points(:,3));

% 求点的距离矩阵
dist_matrix = squareform(pdist(points));

% 求欧式距离矩阵Euclidean distance matrix （注意！欧式距离矩阵一般元素都为距离的平方）
eucl_matrix = dist_matrix.^2;

% 验证其秩是否为小于k+3 = 5
rank(eucl_matrix)

% 添加一行注释测试
rank(eucl_matrix.^0.5)