%% Draw dpm config files

clear;
close all;
clc;

% file name string
fstr = '~/Jamming/CellSim/dpm/pos.test';
bondstr = '~/Jamming/CellSim/dpm/bond.test';

% read in data
dpmData = readMesoPin2D(fstr);

% get number of frames
NFRAMES = dpmData.NFRAMES;

% sim info
NCELLS = dpmData.NCELLS;
nv = dpmData.nv;
L = dpmData.L(1,:);
Lx = L(1);
Ly = L(2);
x = dpmData.x;
y = dpmData.y;
r = dpmData.r;
px = dpmData.px;
py = dpmData.py;
zc = dpmData.zc;
zv = dpmData.zv;
a0 = dpmData.a0;
l0 = dpmData.l0;
t0 = dpmData.t0;
kb = dpmData.kb;

% get preferred shape
calA0 = zeros(NFRAMES,NCELLS);
for ff = 1:NFRAMES
    a0tmp = a0(ff,:);
    l0tmp = l0(ff,:);
    for cc = 1:NCELLS
        p0tmp = sum(l0tmp{cc});
        calA0(ff,cc) = p0tmp^2/(4.0*pi*a0tmp(cc));
    end
end

% particle shape data
p = dpmData.p;
a = dpmData.a;
calA = p.^2./(4.0*pi*a);

% stress data
S = dpmData.S;
P = 0.5*(S(:,1) + S(:,2));

% pulling h
h = dpmData.h;

% draw cell cell contacts
if ~isempty(who('bondstr')) && exist(bondstr,'file')
    % can choose to draw contacts
    drawBonds = 1;

    % load in ctc data
    if drawBonds == 1
        gijList = cell(NFRAMES,1);
        fid = fopen(bondstr);
        for ff = 1:NFRAMES
            NVTOT = sum(nv(ff,:));
            NVVCTCS = 0.5*NVTOT*(NVTOT-1);
            gijtmp = textscan(fid,repmat('%f',1,NVVCTCS),1);
            gijtmp = cell2mat(gijtmp);
            gg = 1;
            gij = zeros(NVTOT);
            for gi = 1:NVTOT
                for gj = (gi+1):NVTOT
                    gij(gi,gj) = gijtmp(gg);
                    gij(gj,gi) = gij(gi,gj);
                    gg = gg + 1;
                end
            end
            gijList{ff} = gij;
        end
        fclose(fid);
    end
else
    drawBonds = 0;
end

% plot area strain
if NFRAMES > 5
    figure(10), clf, hold on, box on;
    dela = a./a0 - 1;
    for nn = 1:NCELLS
        if nn == 1
            plot(h,dela(:,nn),'-','linewidth',2.5,'color',[0.5 0.5 0.5]);
        else
            plot(h,dela(:,nn),'-','linewidth',1.5,'color',[0.5 0.5 0.5]);
        end
    end
    
        p0 = zeros(NFRAMES,NCELLS);
    for ff = 1:NFRAMES
        for nn = 1:NCELLS
            p0(ff,nn) = sum(l0{ff,nn});
        end
    end
    delp = p./p0 - 1;
    for nn = 1:NCELLS
        plot(h,delp(:,nn),'--','linewidth',1.5,'color',[0.5 0.5 0.5]);
    end
    
end
xlabel('$h$','Interpreter','latex');
ylabel('$\Delta_a, \Delta_p$','Interpreter','latex');
ax = gca;
ax.FontSize = 18;

%% Draw cells

% show vertices or not
showverts = 0;

% color by shape or size
colorShape = 4;

if colorShape == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
    cellCLR = jet(NCLR);
elseif colorShape == 2
    % color by preferred shape
    NCLR = 100;
    calA0Bins = linspace(0.999*min(calA0(:)),1.001*max(calA0(:)),NCLR+1);
    cellCLR = jet(NCLR);
elseif colorShape == 3
    clrOrig = [0 0 1];
    clrNew = [1 1 1];
elseif colorShape == 4
    faceColor = 'w';
    NCLR = 100;
    t0All = cell2mat(t0(:));
    mint0 = min(t0All);
    maxt0 = max(t0All);
    t0Bins = linspace(mint0 - 0.001*abs(mint0),maxt0 + 0.001*abs(maxt0),NCLR+1);
    t0clr = jet(NCLR);
    NCBTICKS = 5;
    t0ticks = linspace(0,1,NCBTICKS);
    t0tickLabels = cell(length(t0ticks),1);
    for cc = 1:NCBTICKS
        t0tickLabels{cc} = ['$' sprintf('%0.3g',t0Bins(round(t0ticks(cc)*(NCLR-1)) + 1)) '$'];
    end    
else
    [nvUQ, ~, IC] = unique(nv);
    NUQ = length(nvUQ);
    cellCLROpts = summer(NUQ);
    cellCLR = zeros(NCELLS,3);
    for cc = 1:NCELLS
        cellCLR(cc,:) = cellCLROpts(IC(cc),:);
    end
end

% get frames to plot
if showverts == 0
    FSTART = 1;
    FSTEP = 1;
    FEND = NFRAMES;
%     FEND = FSTART;
else
    FSTART = 1;
    FSTEP = 1;
    FEND = NFRAMES;
%     FEND = FSTART;
end

% make a movie
makeAMovie = 1;
if makeAMovie == 1
    moviestr = 'mesoPin2D_patch2.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
end

fnum = 1;
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    fprintf('printing frame ff = %d/%d\n',ff,FEND);
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    zvf = zv(ff,:);
    nvtmp = nv(ff,:);
    gi = 1;
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        zvtmp = zvf{nn};
        if showverts == 1
            for vv = 1:nvtmp(nn)
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                if zvtmp(vv) < 0
                    rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','c','LineWidth',0.2);
                else
                    rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','b','LineWidth',0.2);
                end
                gi = gi + 1;
            end
        else
            % geometric data
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            vpos = [xtmp, ytmp];
            finfo = [1:nvtmp(nn) 1];
            
            % draw based on color
            if colorShape == 1
                cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
                clr = cellCLR(cbin,:);
                patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
            elseif colorShape == 2
                cbin = calA0(ff,nn) > calA0Bins(1:end-1) & calA0(ff,nn) < calA0Bins(2:end);
                clr = cellCLR(cbin,:);
                patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
            elseif colorShape == 3
                CData = repmat(clrOrig,nvtmp(nn),1);
                CData(zvtmp==-1,:) = repmat(clrNew,sum(zvtmp==-1),1);
                patch('Faces',finfo,'vertices',vpos,CData,'FaceColor','c','EdgeColor','interp');
            elseif colorShape == 4
                t0tmp = t0{ff,nn};
                CData = zeros(nvtmp(nn),3);
                for vv = 1:nvtmp(nn)
                    t0bin = t0tmp(vv) > t0Bins(1:end-1) & t0tmp(vv) < t0Bins(2:end);
                    CData(vv,:) = t0clr(t0bin,:);
                end
                patch('Faces',finfo,'vertices',vpos,'FaceVertexCData',CData,'FaceColor',faceColor,'EdgeColor','interp','linewidth',2.5);
                colormap('jet');
                cb = colorbar;
                cb.Ticks = t0ticks;
                cb.TickLabels = t0tickLabels;
                cb.TickLabelInterpreter = 'latex';
                cb.FontSize = 18;
            else
                patch('Faces',finfo,'vertices',vpos,'FaceColor','b','EdgeColor','k');
            end
        end
    end
    
    % plot box
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [0.2 0.8]*Lx;
    ax.YLim = [0.2 0.8]*Ly;
    
    % plot pin locs
    for cc = 1:NCELLS
        plot(px(ff,cc),py(ff,cc),'ko','markersize',6,'markerfacecolor','k');
    end
    
    % draw bonds if available
    if drawBonds == 1
        NVTOT = sum(nv(ff,:));
        gij = gijList{ff};
        xall = cell2mat(xf');
        yall = cell2mat(yf');
        for gi = 1:NVTOT
            for gj = (gi+1):NVTOT
                if (gij(gi,gj) == 1)
                    xi = xall(gi);
                    yi = yall(gi);
                    
                    xj = xall(gj);
                    yj = yall(gj);
                    
                    plot([xi xj],[yi yj],'k-','linewidth',2);
                end
            end
        end
    end
    
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

%% Draw center particle "growth"

NPLOTS = 5;
NMAX = round(0.75*NFRAMES);
fplots = round(linspace(1,NMAX,NPLOTS));
figure(2), clf, hold on, box on;
for ii = 1:NPLOTS
    ff = fplots(ii);
    xf = x{ff,1};
    yf = y{ff,1};
    rf = r{ff,1};
    zvf = zv{ff,1};
    nvtmp = nv(ff,1);
    
    % get rel pos
    cx = mean(xf);
    cy = mean(yf);
    
    rx = xf - cx;
    ry = yf - cy;
    
    % scale
    ascale = (1 + h(ff))./sqrt(a0(ff,1));
    
    rx = rx.*ascale;
    ry = ry.*ascale;
    rf = rf.*ascale;
    
    % plot
    ip1 = [2:nvtmp 1];
    for vv = 1:nvtmp
        rv = rf(vv);
        xplot = rx(vv) - rv;
        yplot = ry(vv) - rv;
        if zvf(vv) < 0
%             rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clrNew,'LineWidth',0.2);
            plot([rx(vv) rx(ip1(vv))],[ry(vv) ry(ip1(vv))],'-','linewidth',3,'color','b');
        else
%             rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clrOrig,'LineWidth',0.2);
            plot([rx(vv) rx(ip1(vv))],[ry(vv) ry(ip1(vv))],'-','linewidth',1.5,'color','k');
        end
    end
    
    % plot box
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.5 0.5]*Lx;
    ax.YLim = [-0.5 0.5]*Ly;
end