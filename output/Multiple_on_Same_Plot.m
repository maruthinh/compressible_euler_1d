clear all
close all
clc


a1=load('Result.dat');

%a2=load('MoversWEFProb8.dat');
%c3=load('LLFProb6.dat');
%b=load('steadycontact.dat'); %%prob 6
%b1=load('shock_tube.dat');  %%prob1
%b1=load('over_heating.dat');  %%prob2
b1=load('woodward_collela.dat');  %%prob3
%b1=load('shock_collision.dat');  %%prob4
%b1=load('prob5.dat');  %%prob5
%b1=load('prob5.dat');  %%prob5
%b=load('Prob8.dat');  %%prob1
figure = figure('PaperSize',[20.0 20.0],'PaperUnits','inches');
%set(gca,'DefaultTextFontSize',25)
subplot(2,2,1);
plot(a1(:,1), a1(:,2),'k-o','MarkerSize',4,'LineWidth',0.6);
hold on
plot(b1(:,1), b1(:,2),'k','MarkerSize',5,'LineWidth',0.6);
title('Density');
%h = text(0.1,1.1, '(a)');
%set(h, 'FontSize', 15)
%legend('MOVERS-H', 'LLF','Location','northeast')
%legend('boxoff')
%axis([0 1 0.95 1.45]) 

%hold on
subplot(2,2,2)
plot(a1(:,1), a1(:,3),'ko','MarkerSize',4,'LineWidth',0.6)
hold on
plot(b1(:,1), b1(:,3),'k','LineWidth',1.5)
title('Velocity');
%h = text(0.1,0.3, '(b)');
%set(h, 'FontSize', 15)
%legend('MOVERS-H', 'Exact','Location','northeast')
%legend('boxoff')
%axis([0 1 0.05 1.1]) 

subplot(2,2,3)
plot(a1(:,1), a1(:,4),'ko','MarkerSize',4,'LineWidth',0.6)
hold on
plot(b1(:,1), b1(:,4),'k','LineWidth',1.5)
title('Pressure');
%h = text(0.1,3, '(c)');
%set(h, 'FontSize', 15)
%legend('MOVERS-H', 'Exact','Location','northwest')
%legend('boxoff')
%axis([0 1 0 7]) 

subplot(2,2,4)
plot(a1(:,1), a1(:,5),'ko','MarkerSize',4,'LineWidth',0.6)
hold on
plot(b1(:,1), b1(:,5),'k','LineWidth',1.5)
title('Internal energy');
%h = text(0.1,15, '(d)');
%set(h, 'FontSize', 15)
%legend('MOVERS-H', 'Exact','Location','northwest')
%legend('boxoff')
%axis([0 1 0 35]) 

%xlabel('x')
%ylabel('density')
%legend({'MOVERS','Exact'}, 'Position',[0.75,0.85,0.0,0.0])
%legend('MOVERS-LE','Roe','LLF','Exact','Location','northeast')
%legend('MOVERS-LE','LLF','Exact','Location','northeast')
%legend('boxoff')


%set(gca,'FontSize',12)  


%axis([0 1 0.8667 2.7667]) %prob8 
%axis([0 1 0.95 1.45]) %prob6 
%axis([0 1 0.025 1.1])  %prob1
%axis([0 1 0.0 7.0]) %prob5 
%set(gca, 'PaperPosition', [0 0 6 3]);
%axis auto