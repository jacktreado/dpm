%% Draw dpm config files

clear;
close all;
clc;

% file name string
fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
dpmData = readMesoPin2D(fstr);

% get number of frames
NFRAMES = dpmData.NFRAMES;

% sim info
NCELLS = dpmData.NCELLS;
nv = dpmData.nv(1,:);
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
zg = dpmData.zg;
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
elseif colorShape == 2
    % color by preferred shape
    NCLR = 100;
    calA0Bins = linspace(0.999*min(calA0(:)),1.001*max(calA0(:)),NCLR+1);
    cellCLR = jet(NCLR);
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
    FSTART = 10;
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
            vpos = [xtmp, ytmp];
            finfo = [1:nv(nn) 1];
            patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
        end
    end

    % plot pin locs
    for cc = 1:NCELLS
        plot(px(ff,cc),py(ff,cc),'ko','markersize',10,'markerfacecolor','k');
    end
    
    % plot box
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*Lx;
    ax.YLim = [-0.25 1.25]*Ly;
    
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


% print if multiple frames
if NFRAMES > 5
    Sxx = S(:,1);
    Syy = S(:,2);
    Sxy = S(:,3);
    
    figure(10), clf, hold on, box on;
    
    plot(h(Sxx<0),abs(Sxx(Sxx<0)),'ks','markersize',10,'MarkerFaceColor','r');
    plot(h(Sxx>0),Sxx(Sxx>0),'ro','markersize',10);
    
    plot(h(Syy<0),abs(Syy(Syy<0)),'ks','markersize',10,'MarkerFaceColor','b');
    plot(h(Syy>0),Syy(Syy>0),'bo','markersize',10);
    
    xlabel('$h$','Interpreter','latex');
    ylabel('$\Sigma_{xx}$, $\Sigma_{yy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(11), clf, hold on, box on;
    
    plot(h(Sxy<0),abs(Sxy(Sxy<0)),'ks','markersize',10,'MarkerFaceColor','k');
    plot(h(Sxy>0),Sxy(Sxy>0),'ko','markersize',10);
    
    xlabel('$h$','Interpreter','latex');
    ylabel('$\Sigma_{xy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(12), clf, hold on, box on;
    plot(h,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(h,mean(calA,2),std(calA,0,2),'k--','linewidth',2);
    
    xlabel('$h$','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;    
end