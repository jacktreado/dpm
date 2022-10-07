%% Script to draw outcome of adcm2D simulation

clear;
close all;
clc;

% position information
dpmloc = '/Users/jacktreado/Jamming/CellSim/dpm';
posfname = 'adcm2D_test.pos';
posfstr = [dpmloc '/' posfname];

% load into simulation data struct
simdata = readADCM2DTrajectory(posfstr);

% parse data
NFRAMES = simdata.NFRAMES;
NCELLS = simdata.NCELLS;
nv = simdata.nv;
LList = simdata.L;
x = simdata.x;
y = simdata.y;
r = simdata.r;
p = simdata.p;
a = simdata.a;
a0 = simdata.a0;
phi0 = sum(a0,2)./(LList(:,1).*LList(:,2));
phiA = sum(a,2)./(LList(:,1).*LList(:,2));

%% Draw cells over time

% color by shape or size
colorOpt = 0;

% show vertices or not
showverts = 0;

% color option
if colorOpt == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.99,3,NCLR-1);
    calABins = [calABins 10000];
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


% get frames to plot
if showverts == 0
    % single frame
%     FSTART = 1;
%     FEND = FSTART;

    % movie frames
    FSTART = 1;
    FEND = NFRAMES;

    % set step size
    FSTEP = 1;
    DF = FEND - FSTART;
    if DF > 80 && DF <= 400
        FSTEP = 5;
    elseif DF > 400 && DF <= 800
        FSTEP = 10;
    elseif DF > 800
        FSTEP = 20;
    end
else
    FSTART = 10;
    FSTEP = 1;
    FEND = FSTART;
end

% make a movie
makeAMovie = 0;
ctccopy = 0;
if makeAMovie == 1
    moviestr = 'adcm2D_test.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
    ctccopy = -1:1;
end

fnum = 1;
figure(fnum), clf, hold on, box on;
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    fprintf('printing frame ff = %d/%d\n',ff,FEND);
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    L = LList(ff,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        nvtmp = nv(ff,nn);
        
        % get color info
        switch colorOpt
            case 1
                cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
                clr = cellCLR(cbin,:);
            case 2
                cbin = calA0(ff,nn) > calA0Bins(1:end-1) & calA0(ff,nn) < calA0Bins(2:end);
                clr = cellCLR(cbin,:);
            otherwise
                clr = cellCLR(IC(ff,nn),:);
        end
        
        % draw SS polygon
        drawSSPoly(xtmp, ytmp, rtmp, clr);
        
    end
        
    % plot box
    plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*L;
    ax.YLim = [-0.25 1.25]*L;
    
%     colorbar
%     colormap(jet);
%     cb = colorbar;
%     cb.Ticks = [];
    
    
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


