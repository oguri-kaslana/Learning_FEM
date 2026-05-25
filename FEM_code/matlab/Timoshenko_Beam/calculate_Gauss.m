function [Gp, Gw] = calculate_Gauss(number_Gp)
% ============================================================
% calculate_Gauss
%
% 功能：
%   返回一维 Gauss-Legendre 积分点和权重，积分区间 [-1, 1]。
% ============================================================

    if number_Gp == 1

        Gp = 0;
        Gw = 2;

    elseif number_Gp == 2

        Gp = [-sqrt(1/3); sqrt(1/3)];
        Gw = [1; 1];

    elseif number_Gp == 3

        Gp = [-sqrt(3/5); 0; sqrt(3/5)];
        Gw = [5/9; 8/9; 5/9];

    else

        error('calculate_Gauss only supports 1, 2, or 3 Gauss points.');

    end

end
