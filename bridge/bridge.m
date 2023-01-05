clear all; clc; close all;
%%  shear building
tic
load freqtrue_3
load modeltrue_3
load s1_3
load s2_3
%% measurement
n = 15; % problem dimension = No. parameters
y.n = n;
y.freq = freqtrue; y.s1 = s1; % construct measured data y
y.mode = modeltrue; y.s2 = s2;

log_like_fun = 'objfunc';% objective function
%% start TMCMC
A1 = [0.5; ]; B1 = [1.5; ]; % initial values
x_low = [repmat(A1,15,1)];
x_up = [repmat(B1,15,1)]; 

N = 15000;% No. samples
[x,ln_S] = TMCMC_fun(log_like_fun,y,x_low,x_up,N);

toc

%% visualization

% figure (1)
% x=x';
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
npar = 15;
cm = jet(8);
indx=[1 3 6 8 10 11 12 15];
index=1:1:npar;
index(:,indx)=[];
for i=1:8
    set(gca,'colororder',cm)
    plot(x(1:end,indx(i)),'-','LineWidth',1);  
    hold on
end
plot([0,N],[0.9,0.9],'r--','markersize',15,'LineWidth',2);
plot([0,N],[0.8,0.8],'r--','markersize',15,'LineWidth',2);
plot([0,N],[0.7,0.7],'r--','markersize',15,'LineWidth',2);
plot([0,N],[1.2,1.2],'r--','markersize',15,'LineWidth',2);
% legend('\fontsize{15}\bf\theta_1','\fontsize{15}\bf\theta_1_2'...
%         ,'\fontsize{15}\bf\theta_2_5','\fontsize{15}\bf\theta_3_8')
xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5])

figure (4)
for i=1:7
    plot(x(1:end,index(i)),'b-','LineWidth',1);
    hold on
end 
plot([0,N],[1,1],'r--','markersize',15,'LineWidth',2);

xlabel('No. iteration','fontsize',20,'fontname','Times');
ylabel('Stiffness parameters','fontsize',20,'fontname','Times');
set(gca,'fontsize',20);
set(gca,'ylim',[0.5 1.5],'ytick',[0.5:0.2:1.5])
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

MPV = mean(x(0.5*N:end,:));
% std=std(x(0.5*N:end,:));
% save('x')
% set(gca,'xlim',[130 150]);
% set(gca,'xtick',[136.5:0.5:138.5]);
