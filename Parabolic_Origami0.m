clc; clear; close all;

%% 参数设定
a        = 0.005;                   % 抛物线参数 y = a*x^2
diameter = 100;                    % 天线直径
ymax     = a*(diameter/2)^2;       % 抛物面最大高度

%
nY=20;

y        = linspace(0, ymax, nY);               % 水平截面高度 array
nPieces  = 20;                      % 花瓣数目

%% 1. 计算各截面对应的半径
r = sqrt(y./a)

%去掉0，即原点
mask  = (y~=0);         %定义
y2    = y(mask);        %y
r2    = r(mask);        %r
Nlev  = numel(y2);        %层数，y2的维度

%% —— 在这里插入：计算每层的抛物线母线弧长 —— 
% 定义被积函数 (dy/dx = 2*a*x)
arc_fun = @(x) sqrt(1 + (2*a*x).^2);
% 预分配
radiusArc = zeros(Nlev,1);
% 对每一个 r2(k)（对应 x_max = r2(k)）做数值积分
for k = 1:Nlev
    radiusArc(k) = integral( arc_fun, 0, r2(k) );
end

%% 2. 构造 3D 顶点 P(level, piece, xyz)

theta = linspace(0,2*pi,nPieces+1); theta(end)=[];  
P = zeros(Nlev+1, nPieces, 3);
for j = 1:nPieces
    P(1,j,:) = [0,0,0];
end
for k = 1:Nlev
    for j = 1:nPieces
        P(k+1,j,1) = r2(k)*cos(theta(j));
        P(k+1,j,2) = r2(k)*sin(theta(j));
        P(k+1,j,3) = y2(k);
    end
end

%% 3. 计算中央三角形面法向量 n_tri
j0  = 1;
jp  = mod(j0, nPieces)+1;
O   = squeeze(P(1,j0,:))';     % [0 0 0]
A1  = squeeze(P(2,j0,:))';
B1  = squeeze(P(2,jp ,:))';
v1  = A1 - O;
v2  = B1 - O;
n_tri = cross(v1, v2);
n_tri = n_tri / norm(n_tri);

%% 4. 计算中央三角形面与水平面的夹角
% 水平面法向量
n_horiz = [0 0 1];
% 夹角（0°~90°）
angle_tri_horiz = acos( dot(n_tri, n_horiz) ) * 180/pi;

%% 5. 计算各层 trapezoid 面法向量 & 二面角（略，可参照之前代码）

% 5.1. 计算各 facet 的法向量 normals
% --- 中央 triangle facet normal ---
j0 = 1;                        % 任选一个 petal（对称都相同）
jp = mod(j0, nPieces)+1;      % 下一块 petal 的索引
O  = squeeze(P(1,j0,:))';     % [0 0 0]
A1 = squeeze(P(2,j0,:))';     % ring1 左点
B1 = squeeze(P(2,jp ,:))';    % ring1 右点
v1 = A1 - O;
v2 = B1 - O;
n_tri = cross(v1, v2);
n_tri = n_tri / norm(n_tri);

% --- 各层 trapezoid facet normals ---
n_trap = zeros(Nlev-1, 3);
for k = 1:(Nlev-1)
    Pa = squeeze(P(k+1, j0 ,:))';
    Pb = squeeze(P(k+1, jp ,:))';
    Pc = squeeze(P(k+2, j0 ,:))';
    % v_circ: circ direction, v_radial: radial direction
    v_circ   = Pb - Pa;
    v_radial = Pc - Pa;
    nt       = cross(v_circ, v_radial);
    n_trap(k,:) = nt / norm(nt);
end

%5.2. 计算二面角 Dihedral Angles（度数）
% (1) triangle 与第一个 trapezoid
angle_tri_trap = acos( dot(n_tri, n_trap(1,:)) ) * 180/pi;

% (2) 相邻 trapezoid 之间
angle_trap = zeros(Nlev-2,1);
for k = 1:(Nlev-2)
    angle_trap(k) = acos( dot(n_trap(k,:), n_trap(k+1,:)) ) * 180/pi;
end

% 5.3. 显示结果
fprintf('\n【Dihedral Angles】\n');
fprintf(' Triangle ↔ Trap(1): %.3f°\n', angle_tri_trap);
for k = 1:length(angle_trap)
    fprintf(' Trap(%d) ↔ Trap(%d): %.3f°\n', k, k+1, angle_trap(k));
end

%% 6. 输出结果
fprintf('\n【Angle between central triangle and horizontal plane】\n');
fprintf(' Triangle–Horizontal angle = %.3f°\n\n', angle_tri_horiz);

% （可继续输出 triangle vs first trapezoid 及 相邻 trapezoid 间角）
% fprintf(' Triangle ↔ Trap(1): %.3f°\n', angle_tri_trap);
% ...



%% 7. 绘制 3D Parabolic Dish，同时显示每个面的边界线
figure('Color','w'); hold on; axis equal; grid on;
for j = 1:nPieces
    jp = mod(j, nPieces) + 1;
    % ——— 中央三角形 ———
    verts_tri = [...
        squeeze(P(1, j , :))';   % 原点 O
        squeeze(P(2, j , :))';   % 第一环第 j 个点 A1
        squeeze(P(2, jp, :))'];  % 第一环第 jp 个点 B1
    % 使用 patch 并显示边线
    patch( ...
        'Vertices', verts_tri, ...
        'Faces',    [1 2 3],    ...
        'FaceColor','red',      ...
        'FaceAlpha', 0.5,       ...
        'EdgeColor','k',        ... % 黑色边线
        'LineWidth', 1.2        ...
    );
    
    % ——— 各层 trapezoid ———
    for k = 1:(Nlev-1)
        v1 = squeeze(P(k+1, j  , :))';   % 本层左点
        v2 = squeeze(P(k+1, jp , :))';   % 本层右点
        v3 = squeeze(P(k+2, jp , :))';   % 下一层右点
        v4 = squeeze(P(k+2, j  , :))';   % 下一层左点
        
        % 将四个顶点拼成 4×3 数组
        verts_trap = [v1; v2; v3; v4];
        
        % patch 绘制带有边线的四边形
        patch( ...
            'Vertices', verts_trap, ...
            'Faces',    [1 2 3 4], ...
            'FaceColor','green',   ...
            'FaceAlpha', 0.4,      ...
            'EdgeColor','k',       ... % 黑色边线
            'LineWidth', 1.0       ...
        );
        
        % 另外，为了让各条棱特别明显，可以再手动绘制一次边界线
        % （可选，如果觉得 patch 的边线不够显眼）
        % edges = [1 2; 2 3; 3 4; 4 1];
        % for e = 1:4
        %     p1 = verts_trap(edges(e,1), :);
        %     p2 = verts_trap(edges(e,2), :);
        %     plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'k-', 'LineWidth',1.2);
        % end
    end
end

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('3D Parabolic Dish wtih 20 Petals and 20 Radial Segmentations');
view(3);
light; lighting gouraud;

% %% 7. 绘制 3D Parabolic Dish（同之前）
% figure('Color','w'); hold on; axis equal; grid on;
% for j = 1:nPieces
%     jp = mod(j, nPieces) + 1;
%     % 中央三角形
%     verts_tri = [squeeze(P(1,j ,:))'; squeeze(P(2,j ,:))'; squeeze(P(2,jp ,:))'];
%     patch('Vertices', verts_tri, 'Faces', [1 2 3], ...
%           'FaceColor','red','FaceAlpha',0.5,'EdgeColor','k');
%     % 各层 trapezoid
%     for k = 1:(Nlev-1)
%         v1 = squeeze(P(k+1, j   ,:))';
%         v2 = squeeze(P(k+1, jp  ,:))';
%         v3 = squeeze(P(k+2, jp  ,:))';
%         v4 = squeeze(P(k+2, j   ,:))';
%         patch('Vertices',[v1;v2;v3;v4], 'Faces',[1 2 3 4], ...
%               'FaceColor','green','FaceAlpha',0.4,'EdgeColor','none');
%     end
% end
% xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
% title('3D Parabolic Dish with Folding Angles');
% view(3); light; lighting gouraud;





%%





%% —— 新增：计算每一层圆的内接 n 边形边长 —— 
nPolygon = nPieces;  % 你可以改成任意边数，例如 6 边形

side_length = zeros(Nlev,1);
for k = 1:Nlev
    R = r2(k);  % 第 k 层的圆半径
    side_length(k) = 2 * R * sin(pi / nPolygon);
end

% 显示结果
fprintf('\n【每一层圆的内接 %d 边形边长】\n', nPolygon);
for k = 1:Nlev
    fprintf('Layer %d (r=%.3f): Side Length = %.4f\n', k, r2(k), side_length(k));
end


%% 坐标转换

% 假设你已经有了 val 三维数组，大小为 [层数, 列数, 3]
% val(:,:,1) 对应 X， val(:,:,2) 对应 Y， val(:,:,3) 对应 Z

% 提取第一列的所有点 (各层，第1列)
B = squeeze( P(:, 1, :) );   % 维度为 [层数 × 3]，列依次是 [X Y Z]

% 提取第二列的所有点 (各层，第2列)
C = squeeze( P(:, 2, :) );   % 同样是 [层数 × 3]

% 显示结果
disp('B (第一列的 [X Y Z] 坐标)：');
disp(B);
disp('C (第二列的 [X Y Z] 坐标)：');
disp(C);


%% 找梯形角度



for i = 1:nY-1
    b=B(i,:);
    c=C(i,:);
    d=C(i+1,:);
    e=B(i+1,:);


    mid1 = midpoint(c,b);
    
    mid2 = midpoint(d,e);

    if i==1
        mid1m=mid2;
        mid1m(3)=0;
    end
    
    ang(i) = angle(mid2,mid1,mid1m);

    mid1m = mid1;
end




%% 画三角形

% 三个点坐标（3D）%%%%%%%

A = [0, 0, 0]
 b1= B(2,:);
c1 = C(2,:);

allData = [];  % 预先定义一个空矩阵

% 计算边长
AB = norm(b1 - A);
AC = norm(c1 - A);
BC = norm(c1 - b1);

% 打印边长
fprintf('AB = %.4f\n', AB);
fprintf('AC = %.4f\n', AC);
fprintf('BC = %.4f\n', BC);

% 创建 2D 坐标表示：A 在原点，BC 为水平底边，对称分布
% AB 和 AC 是等长的，底边为 BC

% 将 BC 居中放置：B2D 和 C2D 向左右展开
B2D = [-BC/2, sqrt(AB^2 - (BC/2)^2)];
C2D = [ BC/2, sqrt(AC^2 - (BC/2)^2)];
A2D = [0, 0];

% 绘图
figure; hold on; axis equal; grid on;
fill([A2D(1), B2D(1), C2D(1)], [A2D(2), B2D(2), C2D(2)], 'C', 'FaceColor', [0.8, 0.2, 0.2],'FaceAlpha', 0.5);
plot([A2D(1), B2D(1), C2D(1), A2D(1)], [A2D(2), B2D(2), C2D(2), A2D(2)], 'k-o', 'LineWidth', 1.5);

% 标注
text(A2D(1), A2D(2), 'A', 'FontWeight','bold','VerticalAlignment','top');
text(B2D(1), B2D(2), 'b1', 'FontWeight','bold','VerticalAlignment','top');
text(C2D(1), C2D(2), 'c1', 'FontWeight','bold','VerticalAlignment','top');

title('2D 等腰三角形，顶点 A 在原点');
xlabel('X (单位)');
ylabel('Y (单位)');


%% —— 导出 2D 等腰三角形点坐标到 CSV 文件 —— 
triangle_2d = [
    A2D(1), A2D(2), 0;
    B2D(1), B2D(2), 0;
    C2D(1), C2D(2), 0
];

% 写入 CSV（Fusion 360 可导入）
filename = '3角形.csv';
writematrix(triangle_2d, filename);

fprintf('已导出三角形坐标到 CSV 文件：%s\n', filename);



%% 画梯形
f = C2D;


for i = 2:nY-1
    b=B(i,:);
    c=C(i,:);
    d=C(i+1,:);
    e=B(i+1,:);
    
    F = trep(c,b,e,d,f);

    %更新原点位置
    f = F(3,:);

    % 拼成 [X Y 0] 的格式
    csvData = [F, zeros(size(F,1),1)];
    
    % 把这次的结果累加到 allData
    allData = [ allData; csvData ];
end

% 循环结束后一次性写出
writematrix(allData, 'all_trapezoids_flat.csv');
fprintf('已导出所有 trapezoid 点到 CSV：all_trapezoids_flat.csv，共 %d 行\n', size(allData,1));


%%
function [Pts2D] = trep(P1, P2, P3, P4, f)
    % 假设你有四个三维点 P1–P4，且 P1→P2 是短底边，P2 为右下角
    % P1 = [x1, y1, z1];
    % P2 = [x2, y2, z2];
    % P3 = [x3, y3, z3];
    % P4 = [x4, y4, z4];

    % 指定在平面上短底边右下角 P2 映射到的坐标 (u0, v0)
    u0 = f(1);  % 例如 x = 10
    v0 = f(2);  % 例如 y = 5

    %% 1. 构造局部坐标系
    % 1.1 u 轴：沿短底边方向
    u = (P2 - P1) / norm(P2 - P1);

    % 1.2 用 P3 构造平面法向量，再算 v 轴
    n = cross(P2 - P1, P3 - P1);
    n = n / norm(n);
    v = cross(n, u);
    v = v / norm(v);

    %% 2. 将四个点投影到局部 (u,v) 平面并平移
    Pts3D = [P1; P2; P3; P4];
    Pts2D = zeros(4,2);
    for i = 1:4
        vec = Pts3D(i,:) - P1;                % 从 P1 出发的向量
        Pts2D(i,1) = dot(vec, u);             % u 分量
        Pts2D(i,2) = dot(vec, v);             % v 分量
    end

    % 平移，使得 P2 映射到 (u0, v0)
    offset = [u0, v0] - Pts2D(2,:);
    Pts2D = Pts2D + offset;

    %% 3. 在 2D 平面绘制这个梯形（黑色轮廓，绿色填充）
    hold on; 
    axis equal; 
    grid on;

    % 按 P1→P2→P3→P4→P1 顺序构造闭合多边形顶点
    X = Pts2D([1 2 3 4 1], 1);
    Y = Pts2D([1 2 3 4 1], 2);

    % 用 fill 绘制并填充：FaceColor 为绿色，EdgeColor 为黑色
    fill(X, Y,  [0.2, 0.8, 0.2], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);

    % % 如果仍想显示顶点标记，可以再加上 plot（可选）
    % plot(X, Y, 'ks', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    % % 上面 'ks' 表示黑色方块标记，'MarkerFaceColor','g' 让标记内部也是绿色
    % % 若不需要标记，注释掉上面这行即可

    xlabel('X'); 
    ylabel('Y');
    title('2D Plot of a Single Petal of the Unfolded Parabolic Dish');

    hold off;
end



%%
function [M] = midpoint(C, B)
       %% 中点坐标
    
    % % 给定两个三维点
    % C = [10, 20, 5];
    % B = [30, 40, 15];
    
    % 计算中点
    M = (C + B) / 2;
    
    % 显示中点
    fprintf('中点坐标 M = [%.4f, %.4f, %.4f]\n', M);
    
    % % 可视化连线和中点
    % hold on; grid on; axis equal;
    % plot3([C(1), B(1)], [C(2), B(2)], [C(3), B(3)], 'b-', 'LineWidth', 2); % AB线
    % plot3(M(1), M(2), M(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 中点
    % plot3(C(1), C(2), C(3), 'ko', 'MarkerFaceColor', 'k');
    % plot3(B(1), B(2), B(3), 'ko', 'MarkerFaceColor', 'k');
    % 
    % % 标注
    % text(C(1), C(2), C(3), ' C', 'FontSize', 12);
    % text(B(1), B(2), B(3), ' B', 'FontSize', 12);
    % text(M(1), M(2), M(3), ' M (MidPoint)', 'FontSize', 12, 'Color', 'r');
    % 
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    % %title('两点连线与中点 (三维)');
    view(3);

end







%%
function [theta_deg] = angle(M, D, E)
%% 夹角D为顶点
%% 计算向量
v1 = D - M;   % 向量 PD
v2 = E - D;   % 向量 DE

%% 计算夹角（度制）
cos_theta = dot(v1, v2) / (norm(v1) * norm(v2));
theta_deg  = acos( min(max(cos_theta,-1),1) ) * 180/pi;  
% 用 min/max 防止数值误差导致 acos 越界

%% 显示结果
fprintf('向量 PD 与 向量 DE 之间的夹角：%.4f°\n', theta_deg);

%% （可选）三维可视化
% hold on; grid on; axis equal;
% plot3([M(1), D(1)], [M(2), D(2)], [M(3), D(3)], 'b-', 'LineWidth',2);
% plot3([D(1), E(1)], [D(2), E(2)], [D(3), E(3)], 'r-', 'LineWidth',2);
% plot3(M(1), M(2), M(3), 'ko','MarkerFaceColor','k');
% plot3(D(1), D(2), D(3), 'ko','MarkerFaceColor','k');
% plot3(E(1), E(2), E(3), 'ko','MarkerFaceColor','k');
% %text(M(1),M(2),M(3),'  M');
% %text(D(1),D(2),D(3),'  D');
% %text(E(1),E(2),E(3),'  E');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% %title(sprintf('PD 与 DE 之间的夹角 = %.2f°', theta_deg));
view(3);

end