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

%% Draw cells

% show vertices or not
showverts = 0;

% color by shape or size
colorShape = 3;

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
makeAMovie = 0;
if makeAMovie == 1
    moviestr = 'mesoPin2D.mp4';
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
                    rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clrNew,'LineWidth',0.2);
                else
                    rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clrOrig,'LineWidth',0.2);
                end
                gi = gi + 1;
            end
        else
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            vpos = [xtmp, ytmp];
            finfo = [1:nvtmp(nn) 1];
            patch('Faces',finfo,'vertices',vpos,'FaceColor',clrOrig,'EdgeColor','k');
        end
    end
    
    % plot box
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*Lx;
    ax.YLim = [-0.25 1.25]*Ly;
    
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