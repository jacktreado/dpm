%% Draw tumor interface simulations

clear;
close all;
clc;

% file name string
fstr = '~/Jamming/CellSim/dpm/pos.test';
% fstr = 'local/pos/intInit_aN8_ac1.05_tc1.20_aR30_seed130.pos';

% read in data
tumorConfigData = readTumor2DInterface(fstr);

% get number of frames
NFRAMES = tumorConfigData.NFRAMES;

% sim info
NCELLS = tumorConfigData.NCELLS;
tN = tumorConfigData.tN;
nv = tumorConfigData.nv(1,:);
L = tumorConfigData.L;
x = tumorConfigData.x;
y = tumorConfigData.y;
r = tumorConfigData.r;
px = tumorConfigData.px;
py = tumorConfigData.py;
zc = tumorConfigData.zc;
zv = tumorConfigData.zv;
a0 = tumorConfigData.a0;
l0 = tumorConfigData.l0;
t0 = tumorConfigData.t0;
t = tumorConfigData.t;
WP = tumorConfigData.WP;
S = tumorConfigData.S;
P = 0.5*(S(:,1) + S(:,2));

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

% plot wall data
figure(10), plot(t,L(:,1),'k-','linewidth',2);
xlabel('$t$','Interpreter','latex','LineWidth',2);
ylabel('$L_x$','Interpreter','latex','LineWidth',2);
ax = gca;
ax.FontSize = 18;

figure(11), plot(t,WP(:,1),'k-','linewidth',2);
xlabel('$t$','Interpreter','latex','LineWidth',2);
ylabel('$P_{\rm wall}$','Interpreter','latex','LineWidth',2);
ax = gca;
ax.FontSize = 18;

figure(12), plot(t,P,'k-','linewidth',2);
xlabel('$t$','Interpreter','latex','LineWidth',2);
ylabel('$P$','Interpreter','latex','LineWidth',2);
ax = gca;
ax.FontSize = 18;

%% Draw cells

% show vertices or not
showverts = 0;

% color by shape or size
colorShape = 0;

if colorShape == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
    cellCLR = jet(NCLR);
else
    cellCLR = zeros(NCELLS,3);
    cellCLR(1:tN,:) = repmat([1 0 0],tN,1);
    cellCLR((tN+1:end),:) = repmat([1 1 1],NCELLS-tN,1);
end

% get frames to plot
if showverts == 0
    FSTART = 1;
    FSTEP = 10;
    FEND = NFRAMES;
%     FEND = FSTART;
else
    FSTART = NFRAMES;
    FSTEP = 1;
    FEND = FSTART;
end

% make a movie
makeAMovie = 0;
if makeAMovie == 1
    moviestr = 'mesoPin.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
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
    zctmp = zc(ff,:);
    Lx = L(ff,1);
    Ly = L(ff,2);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        clr = cellCLR(nn,:);
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
                for yy = -1:1
                    vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                    finfo = [1:nv(nn) 1];
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
                end
            end
        end
    end


    % plot box
    plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-','linewidth',2);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*L(1,1);
    ax.YLim = [-0.25 1.25]*L(1,2);
    
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