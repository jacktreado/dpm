%% Draw polygons during mesophyll sims

clear;
close all;
clc;

% create file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb01e-3_be50_da0.05_dl7_P1e-4_h0.5_cL0_cB0_seed1.posctc';
% fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.3;
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

% particle shape data
p = mesoData.p(idx,:);
a = mesoData.a(idx,:);
calA = p.^2./(4.0*pi*a);

% stress data
S = mesoData.S(idx,:);
P = 0.5*(S(:,1) + S(:,2));

% construct list of vv contacts
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
polys = cell(NFRAMES,1);
NPOLYS = zeros(NFRAMES,1);
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
    
    % get polygons
    fprintf('* On frame %d, getting void polygons\n',ff);
    cx = cellfun(@mean,x(ff,:))';
    cy = cellfun(@mean,y(ff,:))';
    L = LList(ff,1);
    [mainTiling, ~] = getMesoVoidPolygons(cx,cy,cijtmp,L);
    polys{ff} = mainTiling;
    NPOLYS(ff) = size(polys{ff},1);
    if NPOLYS(ff) >= 6
        fprintf('* On frame %d / %d, got %d void polygons...\n',ff,NFRAMES,NPOLYS(ff));
    else
        fprintf('* On frame %d / %d, got %d void polygons, so ending here...\n',ff,NFRAMES,NPOLYS(ff));
        polys(ff+1:end) = [];
        break;
    end
end
NFRAMES = size(polys,1);


%% Draw cells with void polygons


% show vertices or not
showverts = 0;

% color by shape or size
colorOpt = 1;

if colorOpt == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(1,1.001*max(calA(:)),NCLR+1);
    cellCLR = jet(NCLR);
    cbColorTicks = [0 0.2 0.4 0.6 0.8 1];
    NTICKS = length(cbColorTicks);
    cbTickInds = round(linspace(1,NCLR,NTICKS));
    cbTickLabels = cell(NTICKS,1);
    for tt = 1:NTICKS
        cbTickLabels{tt} = sprintf('$%0.3g$',calABins(cbTickInds(tt)));
    end
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


% get frames to plot
if showverts == 0
    FSTART = 15;
    FSTEP = 1;
%     FEND = NFRAMES-1;
    FEND = FSTART;
else
    FSTART = NFRAMES;
    FSTEP = 1;
    FEND = FSTART;
end

% make a movie
makeAMovie = 0;
ctccopy = 0;
if makeAMovie == 1
%     moviestr = [fpattern '.mp4'];
    moviestr = 'mesoHMin2D_N64_n32_ca1.14_kb01e-3_be50_da0.05_dl7_P1e-4_h0.5_cL0_cB2_seed59.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
    ctccopy = -1:1;
end

fnum = 1;
figure(fnum), clf, hold on, box on;
pclr = jet(7);
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    fprintf('printing frame ff = %d/%d, phi=%0.3g, phi0=%0.3g\n',ff,FEND,phi(ff),phi0(ff));
    
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
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if zctmp(nn) > 0
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',0.2);
                        else
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor',clr,'FaceColor','none');
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
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
                end
            end
        end
    end
    
    if colorOpt == 1
        colormap('jet');
        cb = colorbar;
        cb.Ticks = cbColorTicks;
        cb.TickLabels = cbTickLabels;
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize = 24;
    end
    
    % plot cell-cell contact network
    cijtmp = cijList{ff};
    for ci = 1:NCELLS
        cxi = mean(xf{ci});
        cyi = mean(yf{ci});
        for cj = (ci+1):NCELLS
            if cijtmp(ci,cj) == 1
                cxj = mean(xf{cj});
                cyj = mean(yf{cj});

                dx = cxj - cxi;
                dx = dx - L*round(dx/L);

                dy = cyj - cyi;
                dy = dy - L*round(dy/L);
            
                for xx = -1:1
                    for yy = -1:1
                        plot([cxi cxi + dx] + L*xx,[cyi cyi + dy] + L*yy,'k-','linewidth',2);
                    end
                end
            end
        end
    end
    
%     % plot polygons
%     for pp = 1:length(polys{ff})
%         ptmp = polys{ff}{pp};
%         NE = length(ptmp);
%         if NE < 10
%             patch('Faces',[1:NE 1],'vertices',ptmp,'FaceColor',pclr(NE-2,:),'EdgeColor','w','linewidth',2,'FaceAlpha',1);
%         else
%             patch('Faces',[1:NE 1],'vertices',ptmp,'FaceColor','k','EdgeColor','w','linewidth',2,'FaceAlpha',1);
%         end
%     end
        
    % plot box
    plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*L;
    ax.YLim = [-0.25 1.25]*L;
    
    % if making a movie, save frame
    if makeAMovie == 1
        currframe = getframe(gcf);
        writeVideo(vobj,currframe);
    end
end


% close video object
if makeAMovie == 1
    close(vobj);
end


%% Plot fraction of area contributed by each type of polygon

% polygon area population
minpoly = 3;
maxpoly = 10;
polycheck = minpoly:maxpoly;
NCHECK = length(polycheck);
polyAreaPop = zeros(NFRAMES,NCHECK);
for ff = 1:NFRAMES
    tmppolys = polys{ff};
    NPOLYS = length(tmppolys);
    Abox = LList(ff,1)*LList(ff,2);
    for ii = 1:NPOLYS
        ptmp = tmppolys{ii};
        NE = size(ptmp,1);
        atmp = polyarea(ptmp(:,1),ptmp(:,2));
        if NE < maxpoly
            cc = NE-2;
        else
            cc = maxpoly-2;
        end
        polyAreaPop(ff,cc) = polyAreaPop(ff,cc) + (atmp/Abox);
    end
end


clr = jet(NCHECK);
figure(2), clf, hold on, box on;
for cc = 1:NCHECK
    plot(1-phi(2:end),polyAreaPop(2:end,cc),'-','linewidth',2,'color',clr(cc,:));
end
xlabel('$\varphi - \varphi_{\rm min}$','Interpreter','latex');
ylabel('$f_z(\varphi)$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;
colormap('jet');
cb = colorbar;
cb.Ticks = (polycheck-3)./(NCHECK-1);
cb.TickLabels = {'$3$','$4$','$5$','$6$','$7$','$8$','$9$','$10$','$11$','$\geq 12$'};
cb.TickLabelInterpreter = 'latex';


% also group into 3, 4, 5-7, >7
figure(3), clf, hold on, box on;
clr = jet(4);
for cc = 1:4
    if cc < 3
        plot(1-phi(2:end),polyAreaPop(2:end,cc),'-ko','markerfacecolor',clr(cc,:),'markersize',10);
    elseif cc == 3
        plot(1-phi(2:end),sum(polyAreaPop(2:end,3:5),2),'-ko','markerfacecolor',clr(cc,:),'markersize',10);
    else
        plot(1-phi(2:end),sum(polyAreaPop(2:end,6:end),2),'-ko','markerfacecolor',clr(cc,:),'markersize',10);
    end
end
xlabel('$\varphi - \varphi_{\rm min}$','Interpreter','latex');
ylabel('$f_E$','Interpreter','latex');
ax = gca;
ax.FontSize = 28;
legend({'$E = 3$','$E = 4$','$E = 5-7$','$E \geq 8$'},'Interpreter','latex','FontSize',22,'location','best');