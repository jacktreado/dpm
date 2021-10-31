%% DRAW free meso growth

clear
close all;
clc;

fstr = '~/Jamming/CellSim/dpm/pos.test';
bondstr = '~/Jamming/CellSim/dpm/bond.test';

% read in data
mesoData = readMesoNetwork2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.01;
phi = phi(idx);

% number of frames
NFRAMES = sum(idx);

% sim info
NCELLS = mesoData.NCELLS;
nv = mesoData.nv(idx,:);
L = mesoData.L(1,:);
Lx = L(1);
Ly = L(2);
x = mesoData.x(idx,:);
y = mesoData.y(idx,:);
r = mesoData.r(idx,:);
zc = mesoData.zc(idx,:);
zv = mesoData.zv(idx,:);
a0 = mesoData.a0(idx,:);
l0 = mesoData.l0(idx,:);
t0 = mesoData.t0(idx,:);
kb = mesoData.kb(idx,:);
phi0 = sum(a0,2)/(Lx*Ly);

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
p = mesoData.p(idx,:);
a = mesoData.a(idx,:);
calA = p.^2./(4.0*pi*a);

% stress data
S = mesoData.S(idx,:);
P = 0.5*(S(:,1) + S(:,2));

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

% print if multiple frames
if NFRAMES > 5
    Sxx = S(:,1);
    Syy = S(:,2);
    Sxy = S(:,3);
    
    figure(10), clf, hold on, box on;
    
    plot(phi(Sxx<0),abs(Sxx(Sxx<0)),'ks','markersize',10,'MarkerFaceColor','b');
    plot(phi(Sxx>0),Sxx(Sxx>0),'ro','markersize',10);
    
    plot(phi(Syy<0),abs(Syy(Syy<0)),'ks','markersize',10,'MarkerFaceColor','r');
    plot(phi(Syy>0),Syy(Syy>0),'bo','markersize',10);
    
    xlabel('$1-\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xx}$, $\Sigma_{yy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    figure(11), clf, hold on, box on;
    
    plot(phi(Sxy<0),abs(Sxy(Sxy<0)),'ks','markersize',10,'MarkerFaceColor','k');
    plot(phi(Sxy>0),Sxy(Sxy>0),'ko','markersize',10);
    
    xlabel('$1-\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    figure(12), clf, hold on, box on;
    plot(phi,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(phi,mean(calA,2),std(calA,0,2),'k--','linewidth',2);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    
    
    % compute shear anisotropy
    P = 0.5*(Sxx + Syy);
    sN = (Syy - Sxx)./P;
    sXY = -Sxy./P;
    tau = sqrt(sN.^2 + sXY.^2);
    figure(13), clf, hold on, box on;
    plot(phi,tau,'k-','linewidth',2);
    plot(phi,abs(sN),'b--','linewidth',2);
    plot(phi,abs(sXY),'r--','linewidth',2);
    xlabel('$1-\phi$','Interpreter','latex');
    ylabel('$\hat{\tau}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    
    % plot contact network
    figure(14), clf, hold on, box on;
    plot(phi,zc,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(phi,mean(zc,2),std(zc,0,2),'k--','linewidth',2);
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
end


%% Draw cells

% information for ellipses
th = 0.0:0.01:(2.0*pi);
NTHE = length(th);
ex = cos(th);
ey = sin(th);


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
    t0All = -cell2mat(t0(:));
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


% draw cell cell contacts
if ~isempty(who('ctcstr'))
    % can choose to draw contacts
    drawCTCS = 1;

    % load in ctc data
    if drawCTCS == 1
        cijList = cell(NFRAMES,1);
        zc = zeros(NFRAMES,NCELLS);
        NCTCS = 0.5*NCELLS*(NCELLS-1);
        ltInds =  find(tril(ones(NCELLS),-1));
        frmt = repmat('%f ',1,NCTCS);
        fid = fopen(ctcstr);
        for ff = 1:NFRAMES
            cijtmp = zeros(NCELLS);
            ctmp = textscan(fid,frmt,1);
            ctmp = cell2mat(ctmp);
            cijtmp(ltInds) = ctmp;
            cijtmp = cijtmp + cijtmp';
            cijList{ff} = cijtmp;
            zc(ff,:) = sum(cijtmp > 0,1);
        end
        fclose(fid);
    end
else
    drawCTCS = 0;
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
end

% make a movie
makeAMovie = 0;
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
    fprintf('printing frame ff = %d/%d, phi=%0.3g, phi0=%0.3g\n',ff,FEND,phi(ff),phi0(ff));
    
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
                    rectangle('Position',[xplot, yplot, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','w','LineWidth',0.2);
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
                t0tmp = -t0{ff,nn};
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
    plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-','linewidth',1.5);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.25 1.25]*Lx;
    ax.YLim = [-0.25 1.25]*Ly;
    
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
        currframe = getframe(gca);
        writeVideo(vobj,currframe);
    end
end


% close video object
if makeAMovie == 1
    close(vobj);
end


%% Draw center particle "growth"

NPLOTS = 10;
NMAX = NFRAMES;
fplots = round(linspace(1,NMAX,NPLOTS));
figure(2), clf, hold on, box on;
for ii = 1:NPLOTS
    ff = fplots(ii);
    xf = x{ff,1};
    yf = y{ff,1};
    rf = r{ff,1};
    zvf = zv{ff,1};
    nvtmp = nv(ff,1);
    
    % plot
    ip1 = [2:nvtmp 1];
    for vv = 1:nvtmp
        rv = rf(vv);
        if zvf(vv) <= 0
            plot([xf(vv) xf(ip1(vv))],[yf(vv) yf(ip1(vv))],'-','linewidth',3,'color','b');
        else
            plot([xf(vv) xf(ip1(vv))],[yf(vv) yf(ip1(vv))],'-','linewidth',1.5,'color','k');
        end
    end
    
    % plot box
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
end