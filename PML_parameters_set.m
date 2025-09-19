% --- 参数设置 ---
npml= 2 ;   %PML阶数

PMLwidth = 0.2; % PML layer width in meters

% 高阶 PML (Higher-Order PML) 参数
dx = (max(max(Fscale))*(N)*(N)); % 确保此值与您的网格分辨率一致
dx=0.5*dx;

if(npml==2)
sigma1_opt_val = 0.275 / (150 * pi * dx);
else
    sigma1_opt_val=0;
end
sigma2_opt_val = 2.75   / (150 * pi * dx);
kappa2_opt_val = 0.5;
m1=7;
m2=3;
m3=1;
% 初始化 PML 参数数组
sigma1 = zeros(Np, K);
sigma2 = zeros(Np, K); 
kappa1 = ones(Np, K); 
kappa2 = ones(Np, K); 
alpha1 = zeros(Np, K); 
alpha2 = zeros(Np, K);

% 获取计算域边界
xmin = min(x(:)); xmax = max(x(:));
ymin = min(y(:)); ymax = max(y(:));
zmin = min(z(:)); zmax = max(z(:));

% 为每个节点设置 PML 参数
for k = 1:K
    for i = 1:Np
 
        
        % ========== x-direction PML ==========
        if x(i,k) < xmin + PMLwidth
            xi = (xmin + PMLwidth - x(i,k)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            % 累加来自不同方向的贡献
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib; % 注意：kappa叠加需谨慎
            % alpha2 的总贡献将在所有方向处理完后计算
            
        elseif x(i,k) > xmax - PMLwidth
            xi = (x(i,k) - (xmax - PMLwidth)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib;
        end
        
        % ========== y-direction PML ==========
        if y(i,k) < ymin + PMLwidth
            xi = (ymin + PMLwidth - y(i,k)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib;
            
        elseif y(i,k) > ymax - PMLwidth
            xi = (y(i,k) - (ymax - PMLwidth)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib;
        end
        
        % ========== z-direction PML ==========
        if z(i,k) < zmin + PMLwidth
            xi = (zmin + PMLwidth - z(i,k)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib;
            
        elseif z(i,k) > zmax - PMLwidth
            xi = (z(i,k) - (zmax - PMLwidth)) / PMLwidth;
            sigma1_contrib = sigma1_opt_val * xi^m1;
            sigma2_contrib = sigma2_opt_val * xi^m2;
            kappa2_contrib = kappa2_opt_val * xi^m3;
            
            sigma1(i,k) = sigma1(i,k) + sigma1_contrib;
            sigma2(i,k) = sigma2(i,k) + sigma2_contrib;
            kappa2(i,k) = kappa2(i,k) + kappa2_contrib;
        end
        
        alpha2(i,k) = 0.09 + sigma1(i,k);
        
    end
end

