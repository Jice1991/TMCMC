function [x,ln_S] = TMCMC_fun(log_like_fun,y,x_low,x_up,N)
%
%    log_like_fun: filename for the m-file that computes the log
%                  likelihood.  The m-file must have the format log_like = f(x,y),
%                  where x contains the uncertain variable and y contains
%                  other variables necessary to compute log_like
%    y: other variables necessary to compute the log likelihood
%    x_low & x_up: lower and upper bounds for the uncertain variables
%    x自变量范围
%    N: number of TMCMC samples (per stage)采样个数
%    x: resulting TMCMC samples
%    ln_S: resulting log evidence 
%
%    Note: this TMCMC function assumes that the prior PDF for x is a uniform PDF
%                  defined by x_low & x_up
%
%
%% beta 的取值
b = 0.2; % smaller, flatter orginal:0.2
%% 初始值确定
p = 0; ln_S = 0;%p的初始值为0,
%% drawing prior samples
nnn = length(x_low);%确定中间变量数量
for i=1:nnn,
    x(i,:) = x_low(i) + (x_up(i)-x_low(i))*rand(1,N);%先验样本库
end
for i = 1:N,
    log_like(i,1) = feval(log_like_fun,x(:,i),y);%建立先验样本值
end
%% 计算权重
% 自动选取 p
while p<1,
    low_p = p; up_p = 2; old_p = p;
    while up_p-low_p > 1e-6%确定阈值范围
        current_p = (low_p + up_p)/2;
        temp = exp((current_p-p)*(log_like-max(log_like)));
        cov_temp = std(temp)/mean(temp);%计算协方差
        if cov_temp > 1
            up_p = current_p;
        else
            low_p = current_p;
        end
    end
    p = current_p;%确定p的取值
    if p > 1, break; end % p的终值为1
    weight = temp/sum(temp); % 归一化概率密度比重
    ln_S = ln_S+log(mean(temp))+(p-old_p)*max(log_like);%计算Ln
    %% 计算协方差矩阵
    mu_x = x*weight; sigma = zeros(nnn,nnn);
    for i=1:N,
        sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';%计算协方差矩阵
    end
    sqrt_s = sqrtm(sigma);%求矩阵的平方根
    %% 计算 MCMC
    sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
    for i=1:N,
        now_ind = sam_ind(i);
        x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);%候选样本值
        log_like_c = feval(log_like_fun,x_c,y);%生成函数宏命令
        r = exp(p*(log_like_c - current_log_like(now_ind)));
        if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,%Metropolis准则
            x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
        else
            x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
        end
    end
end
%% 最后一步，也就是p=1时，
temp = exp((1-old_p)*(log_like-max(log_like)));%求临时变量
weight = temp/sum(temp);%求归一化的权重
ln_S = ln_S+log(mean(temp))+(1-old_p)*max(log_like);%中间变量值
mu_x = x*weight; sigma = zeros(nnn,nnn);%期望和方差为转移函数
for i=1:N,
    sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';%求协方差矩阵
end
sqrt_s = sqrtm(sigma);%求协方差矩阵算数平方根
sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
for i=1:N,
    now_ind = sam_ind(i);
    x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);%生成候选值
    log_like_c = feval(log_like_fun,x_c,y);%生成函数宏命令
    r = exp(log_like_c - current_log_like(now_ind));
    if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,%Metropolis准则
        x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;%得到后验分布
    else
        x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
    end
end
