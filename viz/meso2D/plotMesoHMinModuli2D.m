%% Script to plot moduli during simulation
% Note: use only when DM data not included in .hess file

clear;
close all;
clc;


fname = 'mesoDM2D_N32_n32_ca1.14_kb01e-2_be50_da0.02_dl7_P1e-6_h0.5_cL0_cB0_seed100';
fstr = ['local/mesoDM2D_data/' fname '.posctc'];
hessstr = ['local/mesoDM2D_data/' fname '.hess'];

% fstr = '~/Jamming/CellSim/dpm/pos.test';
% hessstr = '~/Jamming/CellSim/dpm/hess.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.375;
phi = phi(idx);

% number of frames
NFRAMES = sum(idx);

% sim info
NCELLS = mesoData.NCELLS;
nv = mesoData.nv(idx,:);
LList = mesoData.L(idx,:);
ctcList = mesoData.ctcs(idx,:);
x = mesoData.x(idx,:);
y = mesoData.y(idx,:);
r = mesoData.r(idx,:);
zc = mesoData.zc(idx,:);
zv = mesoData.zv(idx,:);
a0 = mesoData.a0(idx,:);
l0 = mesoData.l0(idx,:);
t0 = mesoData.t0(idx,:);
kb = mesoData.kb(idx,:);
phi0 = sum(a0,2)./(LList(:,1).*LList(:,2));

% get preferred shape
calA0 = zeros(NFRAMES,NCELLS);
p0 = zeros(NFRAMES,NCELLS);
for ff = 1:NFRAMES
    a0tmp = a0(ff,:);
    l0tmp = l0(ff,:);
    for cc = 1:NCELLS
        p0tmp = sum(l0tmp{cc});
        p0(ff,cc) = p0tmp;
        calA0(ff,cc) = p0tmp^2/(4.0*pi*a0tmp(cc));
    end
end

% particle shape data
p = mesoData.p(idx,:);
a = mesoData.a(idx,:);
calA = p.^2./(4.0*pi*a);

% stress data
S = mesoData.S(idx,:);
P = 0.5*(S(:,1) + S(:,2));


% Read in Hessian data
fid = fopen(hessstr);

% loop over frames, save moduli +  eigenvalues
G = zeros(NFRAMES,1);
B = zeros(NFRAMES,1);
fprintf('Reading in data from Hessian file %s\n',hessstr);
for ff = 1:NFRAMES
    % read in newfr
    newfrstr = fgetl(fid);
    fprintf('%s\n',newfrstr);
    
    % read in packing fraction
    phitmp = textscan(fid,'PACKF %f',1);
    fprintf('PACKF %0.5g\n',phitmp{1});
    
    % read in box size
    Ltmp = textscan(fid,'BOXSZ %f %f',1);
    fprintf('BOXSZ %0.5g\n',Ltmp{1});
    emptystr = fgetl(fid);
    
    % read in shear modulus
    Gtmp = textscan(fid,'SHRMD %f',1);
    fprintf('SHRMD %0.5g\n',Gtmp{1});
    emptystr = fgetl(fid);
    G(ff) = Gtmp{1};
    
    % read in bulk modulus
    Btmp = textscan(fid,'BLKMD %f',1);
    fprintf('BLKMD %0.5g\n',Btmp{1});
    emptystr = fgetl(fid);
    B(ff) = Btmp{1};
    
    % read in newfr
    endfrstr = fgetl(fid);
    fprintf('%s\n',endfrstr);
end
poissonRatio = (B-G)./(B+G);

%% Plot

% color for frames
plotClr = jet(NFRAMES-1);


figure(1), clf, hold on, box on;

yyaxis left;
plot(phi(2) - phi(2:end),B(2:end),'ko','markersize',10,'markerfacecolor','b');
h = ylabel('$B$','Interpreter','latex');
h.Color = 'b';
ax = gca;
ax.FontSize = 22;
ax.YColor = 'b';

yyaxis right;
plot(phi(2) - phi(2:end),G(2:end),'kd','markersize',10,'markerfacecolor','r');
h = ylabel('$G$','Interpreter','latex');
h.Color = 'r';
ax = gca;
ax.FontSize = 22;
ax.YColor = 'r';

xlabel('$\varphi - \varphi_{\rm min}$','Interpreter','latex');




figure(2), clf, hold on, box on;
plot(phi(2) - phi(2:end),poissonRatio(2:end),'-kd','markersize',10,'markerfacecolor','g');
xlabel('$\varphi - \varphi_{\rm min}$','Interpreter','latex');
ylabel('$\nu$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;


% save curves
Gplot = G(2:end);
Bplot = B(2:end);
prPlot = poissonRatio(2:end);
poroPlot = phi(2) - phi(2:end);

kbstr = '1e-2';
bestr = '50';
Pstr = '1e-6';
svstr = ['modCurves_kb' kbstr '_be' bestr '_P' Pstr '.mat'];
save(['local/moduliCurves/' svstr],'Gplot','Bplot','prPlot','poroPlot');
