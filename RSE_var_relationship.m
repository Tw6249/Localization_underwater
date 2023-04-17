% 测试节点数量和通信范围与矩阵补全效果的关系
clear;
clc;

side = 100; % 三维立方体空间的边长

num_max = 30; % 最大节点数量
num_min = 10;
range_max = 100; % 最大通信范围
range_min = 80;
rse_matrix = zeros(num_max-num_min+1, range_max-range_min+1); % 存储归一化重构误差的矩阵

for num_points = num_min:num_max
    points = rand(num_points, 3) * side; % 生成随机点数量
    
    for range = range_min:range_max
        
        dist_matrix = squareform(pdist(points)); % 计算距离矩阵
        dist_matrix_2 = dist_matrix.^2; % 距离矩阵的平方矩阵，一般论文中定义的Euclidean Distance Matrix
        % 不使用随记采样，基于通信距离来选择可获取的距离矩阵元素
        % 通信距离range
        range_2 = range^2; % 通信距离的平方
        
        % 采样元素
        S = dist_matrix_2; 
        S(S<=range_2) = 1;
        S(S>range_2) = 0;
        S = S - diag(diag(S));

        % 计算times次RSE值，然后取平均
        rse_sum = 0;
        times = 5;
        for t = 1:times
            %% 采用不同的方法进行矩阵补全

            % 方法一
            % % 矩阵补全，同时指定补全矩阵为对称矩阵，对角元素为0
            % m = num_points;
            % cvx_begin quiet
            % variable X(m,m)
            % minimize(norm_nuc(X))
            % subject to
            % X.*S==dist_matrix_2.*S;
            % diag(X)==zeros(num_points,1);
            % X == X';
            % cvx_end

            % 方法二：Inexact ALM
            [X,~,~] = inexact_alm_rpca(S);

            


            X = X.^(1/2);
            X = X - diag(diag(X));

            % 归一化重构误差
            rse = norm(X - dist_matrix,'fro')/norm(dist_matrix,'fro');
            rse_sum = rse_sum + rse;
        end
        rse = rse_sum/times;
        rse_matrix(num_points-num_min+1, range-range_min+1) = rse;
    end
end

% 节点密度
rho = (range_min:range_max)/side^3;

% 画图
figure;
[X, Y] = meshgrid(num_min:num_max,range_min:range_max);
surf(X, Y, rse_matrix);
xlabel('节点数量');
ylabel('通信范围');
zlabel('归一化重构误差');
colorbar;

% 画图,密度
figure;
[X, Y] = meshgrid(rho,range_min:range_max);
surf(X, Y, rse_matrix);
xlabel('节点密度（个/m^3）');
ylabel('通信范围（m）');
zlabel('归一化重构误差');
colorbar;

% 等高线
figure;
[X, Y] = meshgrid(rho,range_min:range_max);
contour(X, Y, rse_matrix, [0.05,0.25],'ShowText','on');
xlabel('节点密度（个/m^3）');
ylabel('通信范围（m）');
zlabel('归一化重构误差');
