%% Draw dpm config files

clear;
close all;
clc;

% file name string
% fstr = '~/Jamming/CellSim/dpm/pos.test';
fstr = '~/Jamming/CellSim/dpm/viz/meso2D/local/meso2D_data/meso2D_N64_n24_ca1.06_be5_cL0.01_aL1_cB0_cKb0_seed6.pos';

% read in data
mesoData = readMesoNetwork2D(fstr);

% get number of frames
NFRAMES = mesoData.NFRAMES;

% sim info
NCELLS = mesoData.NCELLS;
nv = mesoData.nv(1,:);
L = mesoData.L(1,:);
Lx = L(1);
Ly = L(2);
x = mesoData.x;
y = mesoData.y;
r = mesoData.r;
zc = mesoData.zc;
zv = mesoData.zv;
zg = mesoData.zg;
a0 = mesoData.a0;
l0 = mesoData.l0;
t0 = mesoData.t0;
kb = mesoData.kb;

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
p = mesoData.p;
a = mesoData.a;
calA = p.^2./(4.0*pi*a);

% stress data
S = mesoData.S;
P = 0.5*(S(:,1) + S(:,2));

% packing fraction
phi = mesoData.phi;

%% Draw cells

% show vertices or not
showverts = 0;

% color by shape or size
colorShape = 1;

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
    cellCLR = winter(NUQ);
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
    FSTART = NFRAMES;
    FSTEP = 1;
    FEND = NFRAMES;
end

% make a movie
makeAMovie = 0;
if makeAMovie == 1
    moviestr = 'mesoNetwork.mp4';
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
        if colorShape == 2
            cbin = calA0(ff,nn) > calA0Bins(1:end-1) & calA0(ff,nn) < calA0Bins(2:end);
            clr = cellCLR(cbin,:);
        elseif colorShape == 1
            cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
            clr = cellCLR(cbin,:);
        else
            clr = cellCLR(IC(nn),:);
        end
        if showverts == 1
            for vv = 1:nv(nn)
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if zctmp(nn) > 0
                            rectangle('Position',[xplot + xx*Lx, yplot + yy*Ly, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',0.2);
                        else
                            rectangle('Position',[xplot + xx*Lx, yplot + yy*Ly, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor',clr,'FaceColor','none');
                        end
                    end
                end
            end
        else
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            for xx = -1:1
                for yy = -1:1
                    vpos = [xtmp + xx*Lx, ytmp + yy*Ly];
                    finfo = [1:nv(nn) 1];
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
                end
            end
        end
    end
    
    % plot box
    plot([0 Lx Lx 0 0], [0 0 Ly Ly 0], 'k-', 'linewidth', 1.5);
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
    
    plot(1-phi(Sxx<0),abs(Sxx(Sxx<0)),'ks','markersize',10,'MarkerFaceColor','r');
    plot(1-phi(Sxx>0),Sxx(Sxx>0),'ro','markersize',10);
    
    plot(1-phi(Syy<0),abs(Syy(Syy<0)),'ks','markersize',10,'MarkerFaceColor','b');
    plot(1-phi(Syy>0),Syy(Syy>0),'bo','markersize',10);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xx}$, $\Sigma_{yy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(11), clf, hold on, box on;
    
    plot(1-phi(Sxy<0),abs(Sxy(Sxy<0)),'ks','markersize',10,'MarkerFaceColor','k');
    plot(1-phi(Sxy>0),Sxy(Sxy>0),'ko','markersize',10);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(12), clf, hold on, box on;
    plot(1-phi,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(1-phi,mean(calA,2),std(calA,0,2),'k--','linewidth',2);
    
    xlabel('$1-\phi$','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;    
end