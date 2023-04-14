% 根据实际情况，依据实际通信范围来确定节点的距离矩阵可获得的元素
% 并且使用距离矩阵的平方矩阵来进行矩阵复原

% 可能在 多维缩放 和 普氏分析 中存在问题

clear;
clc;


% 生成随机点
num_points = 30;
points = rand(num_points, 3) * 100;

% 计算距离矩阵
dist_matrix = squareform(pdist(points));
% 距离矩阵的平方矩阵
dist_matrix_2 = dist_matrix.^2;
% 不使用采样，基于通信距离来选择可行距离元素
% 通信距离为range
range = 82;
range = range^2;
S = dist_matrix_2;
S(S<=range) = 1;
S(S>range) = 0;
S=S-diag(diag(S));

% 矩阵补全，同时指定补全矩阵为对称矩阵，对角元素为0
m = num_points;
cvx_begin
    variable X(m,m)
    minimize(norm_nuc(X))
    subject to 
        X.*S==dist_matrix_2.*S;
        diag(X)==zeros(num_points,1);
        X == X';
cvx_end

X = X.^(1/2);
X = X - diag(diag(X));

% 归一化重构误差
rse = norm(X - dist_matrix,'fro')/norm(dist_matrix,'fro');
rankofd2 = rank(dist_matrix_2);
num_zeros = sum(S(:)==0);

% 将结果转换为坐标
% 多维缩放求其坐标
% points_mds=cmdscale(X,3);
% points_mds=cmdscale(dist_matrix,3)
n=size(X,1);
t=zeros(n,n);
for i=1:n
    for j=1:n
        t(i,j)=-0.5*(X(i,j)^2 -1/n*X(i,:)*X(i,:)' -1/n*X(:,j)'*X(:,j) +1/n^2*sum(sum(X.^2)));
    end
end
[V,D] = eig(t);
points_mds=V(:,1:3)*D(1:3,1:3).^(1/2);

% Procustes Analysis 求绝对坐标
a = 20; % 假设锚点的数量为a
Pa = points(1:a,:); % Pa为锚点的绝对坐标
Pr = points_mds(1:a,:); %Pr为对应锚点行的相对坐标

meanPa = mean(Pa,1); % 中心位置
meanPr = mean(Pr,1);
translation = meanPa' - meanPr'; % 计算位移向量

Pa = Pa - meanPa; % 转移到原点
Pr = Pr - meanPr;

[~,Z,transform] = procrustes(Pa,Pr);
points_mds_abs = points_mds*transform.T;

points_mds_abs = points_mds_abs + ones(n,1)*translation'; % 进行位移
% 普氏分析
% [Pr2,Q] = procrust(Pa,Pr); % 得到旋转矩阵Q 
% points_mds_abs = points_mds*Q; %旋转
% Pr2 = Pr*Q;
% meanPa = mean(Pa,1);
% meanPr = mean(Pr2,1);
% translation = meanPa' - meanPr'; % 计算位移向量
% % Pa = Pa - meanPa; % 转移到原点
% % Pr = Pr - meanPr;
% points_mds_abs = points_mds + ones(n,1)*translation'; % 进行位移


% 绘制图形

% 相对坐标和原始坐标的比较
n = size(points, 1);
n_mds = size(points_mds, 1);

% figure;
% scatter3(points(:,1), points(:,2), points(:,3),'b');
% hold on;
% for i = 1:n
%     text(points(i,1), points(i,2),points(i,3), num2str(i));
% end
% scatter3(points_mds(:,1), points_mds(:,2), points_mds(:,3),'r')
% for i = 1:n
%     text(points_mds(i,1), points_mds(i,2),points_mds(i,3), num2str(i));
% end
% for i = 1:n
%     plot3([points(i,1) points_mds(i,1)], [points(i,2) points_mds(i,2)], [points(i,3) points_mds(i,3)], 'k--');
% end
% hold off;
% legend('points');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Scatter plot of points and points_abs');
% 
% % points和points_mds_abs的行序号
% n = size(points, 1);
% n_mds = size(points_mds_abs, 1);

figure;
scatter3(points(:,1), points(:,2), points(:,3),'b');
hold on;
for i = 1:n
    text(points(i,1), points(i,2),points(i,3), num2str(i));
end
scatter3(points_mds_abs(:,1), points_mds_abs(:,2), points_mds_abs(:,3),'r')
for i = 1:n
    text(points_mds_abs(i,1), points_mds_abs(i,2),points_mds_abs(i,3), num2str(i));
end
for i = 1:n
    plot3([points(i,1) points_mds_abs(i,1)], [points(i,2) points_mds_abs(i,2)], [points(i,3) points_mds_abs(i,3)], 'k--');
end
hold off;
legend('points');
xlabel('x');
ylabel('y');
zlabel('z');
title('Scatter plot of points and points_abs');