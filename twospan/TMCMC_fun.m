function [x,ln_S] = TMCMC_fun(log_like_fun,y,x_low,x_up,N)
%
%    log_like_fun: filename for the m-file that computes the log
%                  likelihood.  The m-file must have the format log_like = f(x,y),
%                  where x contains the uncertain variable and y contains
%                  other variables necessary to compute log_like
%    y: other variables necessary to compute the log likelihood
%    x_low & x_up: lower and upper bounds for the uncertain variables
%    x�Ա�����Χ
%    N: number of TMCMC samples (per stage)��������
%    x: resulting TMCMC samples
%    ln_S: resulting log evidence 
%
%    Note: this TMCMC function assumes that the prior PDF for x is a uniform PDF
%                  defined by x_low & x_up
%
%
%% beta ��ȡֵ
b = 0.2; % smaller, flatter orginal:0.2
%% ��ʼֵȷ��
p = 0; ln_S = 0;%p�ĳ�ʼֵΪ0,
%% drawing prior samples
nnn = length(x_low);%ȷ���м��������
for i=1:nnn,
    x(i,:) = x_low(i) + (x_up(i)-x_low(i))*rand(1,N);%����������
end
for i = 1:N,
    log_like(i,1) = feval(log_like_fun,x(:,i),y);%������������ֵ
end
%% ����Ȩ��
% �Զ�ѡȡ p
while p<1,
    low_p = p; up_p = 2; old_p = p;
    while up_p-low_p > 1e-6%ȷ����ֵ��Χ
        current_p = (low_p + up_p)/2;
        temp = exp((current_p-p)*(log_like-max(log_like)));
        cov_temp = std(temp)/mean(temp);%����Э����
        if cov_temp > 1
            up_p = current_p;
        else
            low_p = current_p;
        end
    end
    p = current_p;%ȷ��p��ȡֵ
    if p > 1, break; end % p����ֵΪ1
    weight = temp/sum(temp); % ��һ�������ܶȱ���
    ln_S = ln_S+log(mean(temp))+(p-old_p)*max(log_like);%����Ln
    %% ����Э�������
    mu_x = x*weight; sigma = zeros(nnn,nnn);
    for i=1:N,
        sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';%����Э�������
    end
    sqrt_s = sqrtm(sigma);%������ƽ����
    %% ���� MCMC
    sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
    for i=1:N,
        now_ind = sam_ind(i);
        x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);%��ѡ����ֵ
        log_like_c = feval(log_like_fun,x_c,y);%���ɺ���������
        r = exp(p*(log_like_c - current_log_like(now_ind)));
        if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,%Metropolis׼��
            x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;
        else
            x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
        end
    end
end
%% ���һ����Ҳ����p=1ʱ��
temp = exp((1-old_p)*(log_like-max(log_like)));%����ʱ����
weight = temp/sum(temp);%���һ����Ȩ��
ln_S = ln_S+log(mean(temp))+(1-old_p)*max(log_like);%�м����ֵ
mu_x = x*weight; sigma = zeros(nnn,nnn);%�����ͷ���Ϊת�ƺ���
for i=1:N,
    sigma = sigma + b^2*weight(i)*(x(:,i)-mu_x)*(x(:,i)-mu_x)';%��Э�������
end
sqrt_s = sqrtm(sigma);%��Э�����������ƽ����
sam_ind = deterministicR((1:N),weight); current_x = x; current_log_like = log_like;
for i=1:N,
    now_ind = sam_ind(i);
    x_c = current_x(:,now_ind) + sqrt_s*randn(nnn,1);%���ɺ�ѡֵ
    log_like_c = feval(log_like_fun,x_c,y);%���ɺ���������
    r = exp(log_like_c - current_log_like(now_ind));
    if r > rand & sum(x_c < x_up)==nnn & sum(x_c > x_low)==nnn,%Metropolis׼��
        x(:,i) = x_c; current_x(:,now_ind) = x_c; current_log_like(now_ind) = log_like_c; log_like(i) = log_like_c;%�õ�����ֲ�
    else
        x(:,i) = current_x(:,now_ind); log_like(i) = current_log_like(now_ind);
    end
end
