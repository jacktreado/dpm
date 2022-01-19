%% Plot frames for model summary figure in eLife paper

clear;
close all;
clc;

% create file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb01e-3_be50_da0.05_dl7_P1e-4_h0.5_cL0_cB2_seed32.posctc';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.35 & phi < 1.0;
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


%% Plot two figures (confluent and network) for comparison

% color by shape or size
colorOpt = 1;

if colorOpt == 1
    % color by real shape
    NCLR = 100;
%     calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
    calABins = linspace(0.99,2.86,NCLR);
    calABins = [calABins 10];
    cellCLR = jet(NCLR);
elseif colorOpt == 2
    % color by preferred shape
    NCLR = 100;
    calA0Bins = linspace(0.999*min(calA0(:)),1.001*max(calA0(:)),NCLR+1);
    cellCLR = jet(NCLR);
else
    [nvUQ, ~, IC] = unique(nv);
    IC = reshape(IC,NFRAMES,NCELLS);
    NUQ = length(nvUQ);
    cellCLR = winter(NUQ);
end

% construct list of contacts
gijList = cell(NFRAMES,1);
for ff = 1:NFRAMES
    nvtmp = sum(nv(ff,:));
    ctctmp = ctcList{ff};
    gijtmp = zeros(nvtmp);
    gi = 1;
    ctchit = 1;
    for ii = 1:nvtmp
        for jj = (ii+1):nvtmp
            if gi == (ctctmp(ctchit)+1)
                gijtmp(ii,jj) = 1;
                gijtmp(jj,ii) = 1;
                ctchit = ctchit + 1;
                if ctchit > length(ctctmp)
                    break;
                end
            end
            gi = gi+1;
        end
        if ctchit > length(ctctmp)
            break;
        end
    end
    gijList{ff} = gijtmp;
end

% get cc contacts
cijList = cell(NFRAMES,1);
for ff = 1:NFRAMES
    gijtmp = gijList{ff};
    cijtmp = zeros(NCELLS);
    nvtmp = nv(ff,:);
    szList = [0 cumsum(nvtmp(1:end-1))] + 1;
    for nn = 1:NCELLS
        for mm = (nn+1):NCELLS
            ctcfound = 0;
            gi = szList(nn);
            for vi = 1:nvtmp(nn)
                gj = szList(mm);
                for vj = 1:nvtmp(mm)
                    if gijtmp(gi,gj) == 1 && ctcfound == 0
                        cijtmp(nn,mm) = 1;
                        cijtmp(mm,nn) = 1;
                        ctcfound = 1;
                    end
                    gj = gj + 1;
                end
                gi = gi + 1;
            end
        end
    end
    cijList{ff} = cijtmp;
end

% list of frames to show
frameList = [1 9];

% show vertices or not
showverts = 0;

% cell to highlight
cell2highlight = 14;

for fff = 1:2
    % reset figure for this frame
    figure(fff), clf, hold on, box on;
    ff = frameList(fff);
    fprintf('printing frame ff = %d, phi=%0.3g, phi0=%0.3g\n',ff,phi(ff),phi0(ff));
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    zctmp = zc(ff,:);
    L = LList(ff,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        cx = mean(xtmp);
        cy = mean(ytmp);
        if nn == cell2highlight
            cxfocus = cx;
            cyfocus = cy;
        end
        rtmp = rf{nn};
        nvtmp = nv(ff,nn);
        if colorOpt == 2
            cbin = calA0(ff,nn) > calA0Bins(1:end-1) & calA0(ff,nn) < calA0Bins(2:end);
            clr = cellCLR(cbin,:);
        elseif colorOpt == 1
            cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
            clr = cellCLR(cbin,:);
        else
            clr = cellCLR(IC(ff,nn),:);
        end
        if showverts == 1
            vpos = [xtmp, ytmp];
            patch('Faces',[1:nvtmp 1],'vertices',vpos,'FaceColor','none','EdgeColor','k','Linewidth',1.5,'LineStyle','-');
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if nn == cell2highlight
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',2);
                        else
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','w','LineWidth',2);
                        end
                    end
                end
            end
        else
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            for xx = -1:1
                for yy = -1:1
                    vpos = [xtmp + xx*L, ytmp + yy*L];
                    finfo = [1:nvtmp 1];
                    if nn == cell2highlight
                        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',2);
                    else
                        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',2);
                    end
%                     text(cx,cy,num2str(nn));
                end
            end
        end
    end
    
    % plot vv contacts
%     gijtmp = gijList{ff};
%     xa = cell2mat(xf');
%     ya = cell2mat(yf');
%     NVTOT = sum(nv(ff,:));
%     for gi = 1:NVTOT
%         xi = xa(gi);
%         yi = ya(gi);
%         for gj = (gi+1):NVTOT
%             if (gijtmp(gi,gj) == 1)
%                 dx = xa(gj) - xi;
%                 dx = dx - L*round(dx/L);
%                 dy = ya(gj) - yi;
%                 dy = dy - L*round(dy/L);
%                 for xx = ctccopy
%                     for yy = ctccopy
%                         plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k--','linewidth',2);
%                     end
%                 end
%             end
%         end
%     end
        
    % plot box
    if showverts == 0
        plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
        axis equal;
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.XLim = [-0.25 1.25]*L;
        ax.YLim = [-0.25 1.25]*L;
    else
        % plot vv contacts
        gijtmp = gijList{ff};
        xa = cell2mat(xf');
        ya = cell2mat(yf');
        NVTOT = sum(nv(ff,:));
        for gi = 1:NVTOT
            xi = xa(gi);
            yi = ya(gi);
            for gj = (gi+1):NVTOT
                if (gijtmp(gi,gj) == 1)
                    dx = xa(gj) - xi;
                    dx = dx - L*round(dx/L);
                    dy = ya(gj) - yi;
                    dy = dy - L*round(dy/L);
                    for xx = -1:1
                        for yy = -1:1
                            plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k-','linewidth',2);
                        end
                    end
                end
            end
        end
    
        axis equal;
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.XLim = [cxfocus - 0.15*L cxfocus + 0.15*L];
        ax.YLim = [cyfocus - 0.15*L cyfocus + 0.15*L];
    end
end

% return;
%% Draw individual cell shape dynamics

% pick random cell
cellidx = 14;

% figure for short-time dynamics
figure(3), clf, hold on, box on;

% get geometric info
xf = x(:,cellidx);
yf = y(:,cellidx);
nvf = nv(:,cellidx);
zvf = zv(:,cellidx);
frList = [2 3 5 6 7 8 9 12];
NFR = length(frList);
LP = LList(end,1);
for ii = 1:NFR
    ff = frList(ii);
    xtmp = xf{ff};
    ytmp = yf{ff};
    zvtmp = zvf{ff};
    nvtmp = nvf(ff);
    
    cx = mean(xtmp);
    cy = mean(ytmp);
    xtmp = xtmp - cx + 0.15*(ii-1)*LP;
    ytmp = ytmp - cy;
    ip1 = [2:nvtmp 1];
    for vv = 1:nvtmp
        if zvtmp(vv) <= 0
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',1.5,'color','r');
        else
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',2.5,'color','k');
        end
    end
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.YLim = [-0.15*LP 0.15*LP];
    ax.XLim = [-0.15*LP NFR*0.15*LP];
end
