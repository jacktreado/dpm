%% just check energy conservation

clear;
close all;
clc;

fstr = '/Users/jacktreado/Jamming/CellSim/dpm/energy.test';
fid = fopen(fstr);
endata = textscan(fid,'%f %f %f %f');
idx = endata{1};
t = endata{2};
K = endata{3};
U = endata{4};

figure(1), clf, hold on, box on;
plot(t,K,'r-','linewidth',2);
plot(t,U,'b-','linewidth',2);
plot(t,U+K,'k-','linewidth',2);
xlabel('time');
ylabel('energy');
ax = gca;
ax.FontSize = 22;
ax.YLim = [0 1.1*max(U+K)];
xl = ax.XLim;
plot(xl,[mean(U+K) mean(U+K)],'k--','linewidth',1.5);