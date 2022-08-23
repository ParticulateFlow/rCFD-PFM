% 2022-08-17 
% Stefan Pirker
% Department of particulate flow modelling
% Linz, Austria

%% clear all

clear all;
close all;
clc;

%% post-proc param's

FS = 12;
MS = 10;
LW = 1;

NHEADERLINES = 0;

%% full CFD

data = importdata('./monitor_rCFD.out',' ',NHEADERLINES);
        
[N_lines, N_columns] = size(data)

figure(1);

plot((data(:,1) - data(1,1)), data(:,2),'r-','LineWidth',LW);

hold on;
grid on;

plot((data(:,1) - data(1,1)), data(:,3),'r:','LineWidth',LW);
plot((data(:,1) - data(1,1)), data(:,4),'r--','LineWidth',LW);

xlabel('Time (s)','FontSize',FS);
ylabel('Mean vertical coord (m)','FontSize',FS);
title('Segregation ','FontSize',FS);
legend('bulk','d = 2 mm','d = 1 mm','Location','northeast');
axis([0 7.0 0.048 0.16]);

%%

figure(2);

hold on;
grid on;

plot((data(:,1) - data(1,1)), data(:,3),'r-','LineWidth',LW);
hold on;
grid on;
plot((data(:,1) - data(1,1)), data(:,4),'r--','LineWidth',LW);

xlabel('Time (s)','FontSize',FS);
ylabel('Mean vertical coord (m)','FontSize',FS);
title('Segregation ','FontSize',FS);
legend('d = 2 mm','d = 1 mm','Location','northeast');
axis([0 7.0 0.048 0.16]);
