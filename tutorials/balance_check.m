#!/usr/bin/octave --silent

close all;
clear;
clc;

% load simulation balance data
casedir = pwd();
filedir = [casedir,'/post/'];
balancefile = [filedir,'balance_monitor.out'];
balancemonitor = csvread(balancefile);

col_t = 1;
col_phase = 2;
col_data = 3;
col_value = 4;
col_target = 5;

ntimesteps = max(balancemonitor(:,col_t))-min(balancemonitor(:,col_t))+1;
stepsize = rows(balancemonitor)/ntimesteps;
nfigures = stepsize;
confidence = 1;

% init figures
hFig(1) = figure;

for ii=1:nfigures

    timesteps = balancemonitor(ii:stepsize:end,col_t);
    iphase = balancemonitor(ii,col_phase);
    idata = balancemonitor(ii,col_data);

    figure(hFig(1));
    clf reset;
    hold on;

    actual_value = balancemonitor(ii:stepsize:end,col_value);
    target_value = balancemonitor(ii:stepsize:end,col_target);

    abs_diff_value = abs(actual_value-target_value);
    epsilon = 0.1*0.5*(max(target_value)-min(target_value));
    delta_k = 0.1*abs_diff_value/epsilon;
    confidence = min(confidence, 1 - max(delta_k));

    plot(timesteps,actual_value,'Color','red','LineStyle','-');
    plot(timesteps,target_value,'Color','black','LineStyle','-');
    plot(timesteps,target_value+epsilon,'Color','black','LineStyle','--');
    plot(timesteps,target_value-epsilon,'Color','black','LineStyle','--');

    xlim auto;
    ylim auto;
    xlabel('time step');
    grid on;

    print(hFig(1),[filedir,'balance_phase',int2str(iphase),'_data',int2str(idata),'.eps'],'-color','-deps');
end

printf("%i\n", floor(100*confidence+0.5));

