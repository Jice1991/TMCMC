clear all; clc; close all;
%%  shear building
tic
load freqtrue_exp8
load modeltrue_exp8
load s1_exp8
load s2_exp8
%% measurement
n = 40; % problem dimension = No. parameters
y.n = n;
y.freq = freqtrue; y.s1 = s1; % construct measured data y
y.mode = modeltrue; y.s2 = s2;

log_like_fun = 'objfunc';% objective function
%% start TMCMC
A1 = [0.7; 0.7]; B1 = [1.2; 1.2]; % initial values
x_low = [repmat(A1,20,1)];
x_up = [repmat(B1,20,1)]; 

N = 40000;% No. samples
[x,ln_S] = TMCMC_fun(log_like_fun,y,x_low,x_up,N);

toc

%% visualization

% figure (1)
x=x';
% cm = jet(n);
% for i=1:n
%     set(gca,'colororder',cm)
%     plot(x(1:end,i),'LineWidth',1.2); 
%     hold on
% end
% ylim([0.7 1.2]);
% set(gca,'ytick',[0.7:0.1:1.2]);
% set(gca,'fontsize',20);
% 
% figure(2) 
% for i=1:n
%     subplot(10,4,i)
%     h=histfit(x(0.3*N:end,i)); %burning 30%
%     h(1).FaceColor = [0 0 1]; 
%     h(2).Color = [1 0 0];
%     set(gca,'fontsize',20);
%     hold on
% end

figure(3)
plot(x(1:end,1),'m-','LineWidth',1); 
    hold on
plot(x(1:end,12),'g-','LineWidth',1); 
    hold on
plot(x(1:end,25),'y-','LineWidth',1); 
    hold on
plot(x(1:end,38),'k-','LineWidth',1); 
hold on
% plot([0,n],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
% plot([0,n],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_1_2'...
%         ,'\fontsize{15}\bf\theta_2_5','\fontsize{15}\bf\theta_3_8')
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.7 1.2],'ytick',[0.7:0.1:1.2])

figure (4)
for i=1:40
    if i==1
        plot(x(1:end,i),'m-','LineWidth',1); 
    hold on
    elseif i==12
            plot(x(1:end,i),'g-','LineWidth',1); 
    hold on
    elseif i==25
            plot(x(1:end,i),'y-','LineWidth',1); 
    hold on
    elseif i==38
            plot(x(1:end,i),'k-','LineWidth',1); 
    hold on
    else 
        plot(x(1:end,i),'b-','LineWidth',1); 
    end
%     plot(n,y(i),'x','markersize',10,'LineWidth',2);
end
hold on
% plot([0,n],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
% plot([0,n],[0.8,0.8],'r--','markersize',15,'LineWidth',2);

xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.7 1.2],'ytick',[0.7:0.1:1.2])
% legend('Samples','Gaussian PDF','Location','Northeast');
% legend('boxon')
% ylabel('Number of samples');
% xlabel('x_1','fontsize',18,'fontweight','bold');


% subplot(1,3,2)
% histfit(x(0.3*N:end,2)) %burning 30%
% 
% xlabel('x_2','fontsize',18,'fontweight','bold');
% ylabel('Number of samples');
% set(gca,'fontsize',20);
% 
% subplot(1,3,3)
% histfit(x(0.3*N:end,3)) %burning 30%
% 
% xlabel('x_3','fontsize',18,'fontweight','bold');
% ylabel('Number of samples');
% set(gca,'fontsize',20);

MPV = mean(x(0.5*N:end,1:n));
std=std(x(0.5*N:end,1:n));
save('x')
% set(gca,'xlim',[130 150]);
% set(gca,'xtick',[136.5:0.5:138.5]);
