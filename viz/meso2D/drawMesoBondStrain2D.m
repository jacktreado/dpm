%% Script to plot fraction of compressed vs extended bonds in given simulation

clear;
close all;
clc;

% create file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb0.4_be100_h1_da0.5_dl5_cL0.5_cB4_t0m0.2_P1e-3_seed14.posctc';
% fstr = 'local/mesoDM2D_data/mesoDM2D_N32_n32_ca1.14_kl1_kb01e-3_be50_da0.02_dl10_P1e-4_seed27.posctc';
% fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.45;
NSKIP = 0;
idx(1:NSKIP) = zeros(NSKIP,1);
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
a = mesoData.a(idx,:);
a0 = mesoData.a0(idx,:);
l0 = mesoData.l0(idx,:);
t0 = mesoData.t0(idx,:);
kb = mesoData.kb(idx,:);
phi0 = sum(a0,2)./(LList(:,1).*LList(:,2));
phiA = sum(a,2)./(LList(:,1).*LList(:,2));

% loop over frames, get patch phi
phipatch = zeros(NFRAMES,1);
for ff = 1:NFRAMES
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    apatch = zeros(NCELLS,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        cx = mean(xtmp);
        cy = mean(ytmp);
        rx = xtmp - cx;
        ry = ytmp - cy;
        rads = sqrt(rx.^2 + ry.^2);
        xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
        ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
        apatch(nn) = polyarea(xtmp,ytmp);
    end
    phipatch(ff) = sum(apatch)/(LList(ff,1)*LList(ff,2));
end
dporo = max(phipatch) - phipatch;
FCMP = NFRAMES;

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
P = mesoData.P(idx);
Sxx = mesoData.S(idx,1);
Syy = mesoData.S(idx,2);
Pvirial = 0.5*(Sxx + Syy);

% angle data
t0mean = cellfun(@mean,t0);
t0min = cellfun(@min,t0);
t0max = cellfun(@max,t0);




% compute fraction of compressed vs strained bonds
fstretch = zeros(NFRAMES,1);
for ff = 1:NFRAMES
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    gijtmp = gijList{ff};
    xa = cell2mat(xf');
    ya = cell2mat(yf');
    ra = cell2mat(rf');
    NVTOT = sum(nv(ff,:));
    nstretch = 0;
    for gi = 1:NVTOT
        xi = xa(gi);
        yi = ya(gi);
        ri = ra(gi);
        for gj = (gi+1):NVTOT
            dx = xa(gj) - xi;
            dx = dx - L*round(dx/L);
            dy = ya(gj) - yi;
            dy = dy - L*round(dy/L);
            bl = sqrt(dx*dx + dy*dy);
            sij = ri + ra(gj);
            if (gijtmp(gi,gj) == 1)
                for xx = ctccopy
                    for yy = ctccopy
                        plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'r-','linewidth',2.5);
                    end
                end
                nstretch = nstretch + 1;
            elseif bl < sij
                for xx = ctccopy
                    for yy = ctccopy
                        plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k-','linewidth',2.5);
                    end
                end
            end
        end
    end
end













