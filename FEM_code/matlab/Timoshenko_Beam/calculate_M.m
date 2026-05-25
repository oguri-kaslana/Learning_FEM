function M = calculate_M(E, I, kappa)
% ============================================================
% calculate_M
%
% 功能：
%   根据曲率计算弯矩。
%
% 说明：
%   公式仍为 M = E*I*kappa。
%   与 Euler 梁不同的是，Timoshenko 梁中的 kappa = dtheta/dx，
%   而不是 d²w/dx²。
% ============================================================

    M = E * I * kappa;

end
