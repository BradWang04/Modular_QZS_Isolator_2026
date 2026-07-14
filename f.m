% 定义二阶微分方程的函数 f
function res = f(t, y, v)
    % 这里定义你的二阶微分方程
    res = -10 * y - 2 * v; % 示例方程 y'' = -10y - 2v
end