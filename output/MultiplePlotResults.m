clear all
close all
clc

a1=load('MoversLEProb6.dat');
%a2=load('MoversWEFProb8.dat');
a3=load('LLFProb6.dat');
%b=load('steadycontact.dat'); %%prob 6
%b=load('shock_tube.dat');  %%prob1
%b=load('woodward_collela.dat');  %%prob3
%b=load('shock_collision.dat');  %%prob4
b=load('MoversLEProb6.dat');  %%prob5
%b=load('Prob8.dat');  %%prob1
figure = figure('PaperSize',[20.0 20.0],'PaperUnits','inches');
%set(gca,'DefaultTextFontSize',25)
plot(a1(:,1), a1(:,2),'k*','MarkerSize',8,'LineWidth',0.6)
hold on
%plot(a2(:,1), a2(:,2),'k+','MarkerSize',5,'LineWidth',0.6)
%hold on
plot(a3(:,1), a3(:,2),'ko','MarkerSize',4,'LineWidth',0.6)
hold on
plot(b(:,1), b(:,2),'k','LineWidth',1.5)
xlabel('x')
ylabel('density')
%legend({'MOVERS','Exact'}, 'Position',[0.75,0.85,0.0,0.0])
%legend('MOVERS-LE','Roe','LLF','Exact','Location','northeast')
legend('MOVERS-LE','LLF','Exact','Location','northeast')
legend('boxoff')

h = text(0.1,1, '(a)');
set(h, 'FontSize', 20)
%set(gca,'FontSize',12)  


%axis([0 1 0.8667 2.7667]) %prob8 
axis([0 1 0.95 1.45]) %prob6 
%axis([0 1 0.025 1.1])  %prob1
%axis([0 1 0.0 7.0]) %prob5 
%set(gca, 'PaperPosition', [0 0 6 3]);
%axis auto