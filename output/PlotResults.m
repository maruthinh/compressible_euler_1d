clear all
close all
clc

a=load('Result.dat');
%b=load('steadycontact.dat'); %%prob 6
%b=load('shock_tube.dat');  %%prob1
%b=load('woodward_collela.dat');  %%prob3
%b=load('shock_collision.dat');  %%prob4
%b=load('prob5.dat');  %%prob5
%b=load('Prob8.dat');  %%prob1
%b=load('over_heating.dat');  %%prob2
%b=load('over_heating.dat');  %%prob2
%b=load('prob9.out');  %%prob9
%b=load('prob9_at_t_1.out');
%b=load('rand_choice_prob10_026.out');%%prob10 at t=0.026    
%b=load('rand_choice_prob10_038.out'); %%prob10 at t=0.038
b=load('prob11_ref_10000.dat'); %%prob12 at t=1.8
%b=load('prob12_ref_50000.dat'); %%prob12 at t=1.8
%c=load('prob9_llf.dat'); %%prob 6
figure = figure('PaperSize',[20.0 20.0],'PaperUnits','inches');
set(gca,'DefaultTextFontSize',25)
plot(a(:,1), a(:,2),'ko','MarkerSize',3,'LineWidth',0.5)
hold on
plot(b(:,1), b(:,2),'k','LineWidth',1.0)
%hold on 
%plot(c(:,1), c(:,2),'ko','LineWidth',0.8)
xlabel('x')
ylabel('density')
%legend({'MOVERS','Exact'}, 'Position',[0.75,0.85,0.0,0.0])
legend('MOVERS-H','Reference', 'LLF', 'Location','northwest')
legend('boxoff')

%h = text(0.9,1.4, '(a)');
%set(h, 'FontSize', 20)
%set(gca,'FontSize',12)  

%axis([0 1 0.8667 2.7667]) %prob8 
%axis([0 1 0.95 1.45]) %prob6 
%axis([0 1 0.025 1.1])  %prob1
%axis([0 1 0.0 6.5])  %prob3
%axis([0 1 2.0 35.0])  %prob4
%axis([0 1 0.0 6.5])  %prob4
%axis([0.4 0.9 -0.2 6.7])  %prob4
%set(gca, 'PaperPosition', [0 0 6 3]);
%axis auto





