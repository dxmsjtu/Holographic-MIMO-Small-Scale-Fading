%This Matlab script can be used to generate the variances of the Fourier
%random coefficients in the Fourier plane-wave series expansion of a 2D
%channel in Eq.(43). The script is valid for isotropic channels only and is 
%based on the theoretical computation in Appendix IV.C (part I) of the article:
%
%A. Pizzo, T. L. Marzetta and L. Sanguinetti, "Spatially-Stationary Model
%for Holographic MIMO Small-Scale Fading," in IEEE Journal on Selected Areas
%in Communications, vol. 38, no. 9, pp. 1964-1979, Sept. 2020,
%doi: 10.1109/JSAC.2020.3000877.
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
clear all; close all;   %clc;
%% Parameters
%array size in number of wavelenghts (must be integer)
Lx = 16;

%% Variances of Fourier random coefficients
%discrete wavenumber frequencies
l_vec = [-Lx:1:Lx-1]';

%compute Fourier variances (2*Lx vector)
variances = asin((l_vec+1)/Lx) - asin(l_vec/Lx);

%normalize variances to their maximum
variances_norm = variances/max(variances);

%plot the variances in dB within the support segment
figure;FontSize =28;
plot(l_vec,10*log10(variances_norm));
xlabel('$\ell$','Interpreter','Latex');
ylabel('$\sigma^2_{\ell}$ (dB)','Interpreter','Latex');
xlim([-Lx Lx])
grid on; box on;
set(gca,'FontSize',FontSize);
set(gcf, 'Position', get(0, 'Screensize'));
% % 生成随机信号
% t = 0:0.1:10; % 时间向量
% x = sin(2*pi*0.5*t) + 0.5*randn(size(t)); % 添加噪声的正弦信号
% 
% % 计算相关矩阵
% C = cov(x); % 计算协方差矩阵
% 
% % 计算特征值和特征向量
% [V, D] = eig(C); % V 中包含特征向量，D 中包含特征值
% 
% % 选择主要成分
% coeff = V' * x; % 计算主要成分系数
% 
% % 重构信号
% x_reconstructed = V * coeff;
% 
% % 绘制原始信号和重构信号
% plot(t, x, 'b', t, x_reconstructed, 'r--');
% legend('原始信号', '重构信号');
% xlabel('时间');
% ylabel('幅值');
% title('K-L级数重构示例');

