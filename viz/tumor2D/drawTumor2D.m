%% Draw tumor interface simulations

clear;
close all;
clc;

% file name string
simtype     = 'mono';
NTstr       = '2e6';
Nstr        = '64';
nstr        = '24';
tcstr       = '1.2';
l1str       = '0.01';
Drstr       = '0.1';
gamttstr    = '0';
taustr      = '1';
seedstr     = '12';

fpattern = [simtype '_NT' NTstr '_N' Nstr '_n' nstr '_tc' tcstr '_l1' l1str '_Dr' Drstr '_gamtt' gamttstr '_tau' taustr '_seed' seedstr];
floc = 'local/pos';
fstr = [floc '/' fpattern '.pos'];
fstr = '~/Jamming/CellSim/dpm/pos.test';

% use to never run movie first
initClear = 1;

% read in data
tumorConfigData = readTumor2D(fstr);

% get number of frames
NFRAMES = tumorConfigData.NFRAMES;

% sim info
NCELLS = tumorConfigData.NCELLS;
tN = tumorConfigData.tN;
nv = tumorConfigData.nv(1,:);
L = tumorConfigData.L(1,:);
Lx = L(1);
Ly = L(2);
x = tumorConfigData.x;
y = tumorConfigData.y;
r = tumorConfigData.r;
zc = tumorConfigData.zc;
zv = tumorConfigData.zv;
a0 = tumorConfigData.a0;
l0 = tumorConfigData.l0;
t0 = tumorConfigData.t0;
psi = tumorConfigData.psi;

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
p = tumorConfigData.p;
a = tumorConfigData.a;
calA = p.^2./(4.0*pi*a);

%% Draw cells

% show vertices or not
showverts = 0;

% make a movie
makeAMovie = 0;

% color by (0=red, 1=shape, 2=cos(psi))
colorOpt = 0;

% color by real shape
if colorOpt == 1
    NCLR = 100;
    calABins = linspace(1,1.1*max(calA(:)),NCLR+1);
    cellCLR = jet(NCLR);
elseif colorOpt == 2
    NCLR = 100;
    cpsiBins = linspace(-1,1,NCLR+1);
    cellCLR = winter(NCLR);
else
    [nvUQ, ~, IC] = unique(nv);
    NUQ = length(nvUQ);
    cellCLR = jet(4*NUQ);
end

% get frames to plot
if showverts == 0
    if initClear == 1
        FSTART = NFRAMES;
        FSTEP = 1;
        FEND = FSTART;
    else
        FSTART = 1;
        FSTEP = 1;
        FEND = NFRAMES;
    end
else
    FSTART = 4;
    FSTEP = 1;
    FEND = FSTART;
end

% check whether or not to actually make movie
if initClear == 1
    makeAMovie = 0;
    initClear = 0;
end

% make a movie
if makeAMovie == 1
    moviestr = [fpattern '.mp4'];
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
end

fnum = 1;
figure(fnum), clf, hold on, box on;
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    f = gcf;
    f.Color = 'w';
    fprintf('printing frame ff = %d/%d\n',ff,FEND);
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    zctmp = zc(ff,:);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        if colorOpt == 1
            calAtmp = calA(ff,nn);
            calAInd = calAtmp < calABins(2:end) & calAtmp > calABins(1:end-1);
            clr = cellCLR(calAInd,:);
        elseif colorOpt == 2
            cpsitmp = cos(psi(ff,nn));
            cpsiind = cpsitmp < cpsiBins(2:end) & cpsitmp > cpsiBins(1:end-1);
            clr = cellCLR(cpsiind,:);
        else
            clr = cellCLR(IC(nn),:);
        end
        if showverts == 1
            for vv = 1:nv(nn)
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',0.2);
            end
        else
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            for xx = 0
                for yy = 0
                    vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                    finfo = [1:nv(nn) 1];
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
                end
            end
        end
    end


    % plot box
%     plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-','linewidth',2);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*Lx;
    ax.YLim = [-0.25 1.25]*Ly;
    
    % if making a movie, save frame
    if makeAMovie == 1
        currframe = getframe(f);
        writeVideo(vobj,currframe);
    end
end


% close video object
if makeAMovie == 1
    close(vobj);
end