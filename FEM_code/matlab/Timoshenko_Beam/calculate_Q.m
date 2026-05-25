function Q = calculate_Q(kappa_s, G, A, gamma)
% ============================================================
% calculate_Q
%
% 功能：
%   根据剪切应变 gamma 计算 Timoshenko 梁剪力 Q。
%
% 说明：
%   Euler-Bernoulli 梁中通常不直接由剪切应变计算剪力；
%   Timoshenko 梁中 Q = kappa_s * G * A * gamma。
% ============================================================

    Q = kappa_s * G * A * gamma;

end
