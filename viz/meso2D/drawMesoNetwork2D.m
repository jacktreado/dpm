%% Draw dpm config files
% NOTE: need to get configs with contact info

clear;
close all;
clc;

% create file name

% parameters
Nstr = '32';
nstr = '32';
castr = '1.08';
kb0str = '1e-3';
bestr = '6';
hstr = '0.2';
cLstr = '5';
aLstr = '1';
cBstr = '0';
cKbstr = '0';

% seed
seed = 5;
seedstr = num2str(seed);

% file name str
floc = '~/Jamming/CellSim/dpm/viz/meso2D/local/meso2D_data';
% fpattern = ['meso2D_N' Nstr '_n' nstr '_ca' castr '_be' bestr '_cL' cLstr '_aL' aLstr '_cB' cBstr '_cKb' cKbstr '_seed' seedstr];
% fpattern = ['meso2D_N' Nstr '_n' nstr '_ca' castr '_kb0' kb0str '_be' bestr '_cL' cLstr '_aL' aLstr '_cB' cBstr '_seed' seedstr];
fpattern = ['meso2D_N' Nstr '_n' nstr '_ca' castr '_kb0' kb0str '_be' bestr '_h' hstr '_cL' cLstr '_aL' aLstr '_cB' cBstr '_cKb' cKbstr '_seed' seedstr];
fstr = [floc '/' fpattern '.pos'];
fstr = '~/Jamming/CellSim/dpm/pos.test';

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
LList = mesoData.L;
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

% stress data
S = mesoData.S(idx,:);
P = 0.5*(S(:,1) + S(:,2));

% print if multiple frames
if NFRAMES > 5
    Sxx = S(:,1);
    Syy = S(:,2);
    Sxy = S(:,3);
    P = 0.5*(Sxx + Syy);
    
    figure(10), clf, hold on, box on;
    
    plot(1-phi,P,'ko','markersize',10,'markerfacecolor','b');
%     plot(1-phi,P,'-b','linewidth',2);
    xlabel('$1-\phi$','Interpreter','latex');
    ylabel('$P$','Interpreter','latex');
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
    plot(1:NFRAMES,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(calA,2),'k-','linewidth',3);
    xlabel('frame id','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    
    
    % compute shear anisotropy
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
    plot(1:NFRAMES,zc,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(zc,2),'k--','linewidth',2.5);
    xlabel('frame id.','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    % plot area deviations
    ea = 0.5*(a./a0 - 1).^2;
    figure(15), clf, hold on, box on;
    plot(phi,ea,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(phi,mean(ea,2),std(ea,0,2),'k--','linewidth',2);
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$U_a$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    ax.YScale = 'log';
    
     % plot packing fractions
    figure(16), clf, hold on, box on;
    yyaxis left
    plot(1:NFRAMES,phi,'ko','markersize',10);
    ylabel('$\phi$','Interpreter','latex');
    yyaxis right
    plot(1:NFRAMES,P,'ks','markersize',10);
    ylabel('$P$','Interpreter','latex');
    xlabel('frame','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
     % plot p0 vs phi0
    figure(17), clf, hold on, box on;
    plot(mean(calA0,2),phi,'ko','markersize',10);
    xlabel('$\mathcal{A}_0$','Interpreter','latex');
    ylabel('$\phi$','Interpreter','latex');
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
showverts = 1;

% color by shape or size
colorOpt = 1;

if colorOpt == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
    calA_bc = 0.5*(calABins(1:end-1) + calABins(2:end));
    cellCLR = jet(NCLR);
    NCLRBAR = 5;
    CLRINDS = 1:round(NCLR/NCLRBAR):NCLR;
    cbTickCell = cell(NCLRBAR+1,1);
    cbTickCell{1} = '$1$';
    for cc = 2:NCLRBAR
        calAtmp = calA_bc(CLRINDS(cc));
        cbTickCell{cc} = ['$' sprintf('%0.3g',calAtmp) '$'];
    end
    cbTickCell{end} = ['$' sprintf('%0.3g',calA_bc(end)) '$'];
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
    FEND = FSTART;
end

% make a movie
makeAMovie = 1;
if makeAMovie == 1
    moviestr = [fpattern '.mp4'];
%     moviestr = 'meso2D_enthalpy.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
end

fnum = 1;
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    f = figure(fnum);
    f.Color = 'w';
    clf, hold on, box on;
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
        
%         % draw effective radius
%         rx = xtmp - cx;
%         ry = ytmp - cy;
%         rads = sqrt(rx.^2 + ry.^2);
%         xedge = xtmp + 1*rtmp.*(rx./rads);
%         yedge = ytmp + 1*rtmp.*(ry./rads);
%         ex = xedge - cx;
%         ey = yedge - cy;
%         r2 = sum(ex.^2 + ey.^2)/nvtmp;
%         reff = sqrt(r2);
%         for xx = -1:1
%             for yy = -1:1
%                 rectangle('Position',[cx - reff + xx*L, cy - reff + yy*L, 2*reff, 2*reff],'Curvature',[1 1],'FaceColor',clr);
%             end
%         end
        
        if showverts == 1
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'linewidth',2);
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
        
%         text(cx,cy,num2str(nn));
%         plot(mean(xtmp),mean(ytmp),'wo','markersize',8,'markerfacecolor','w');
        
%         % get shape tensor
%         cx = mean(xtmp);
%         cy = mean(ytmp);
%         rx = xtmp - cx;
%         ry = ytmp - cy;
%         rn = sqrt(rx.^2 + ry.^2);
%         urx = rx./rn;
%         ury = ry./rn;
%         Gxx = sum(rx.*urx)/nvtmp;
%         Gyy = sum(ry.*ury)/nvtmp;
%         Gxy = sum(rx.*ury)/nvtmp;
%         [V,D] = eig([Gxx, Gxy; Gxy, Gyy]);
%         lambda = diag(D);
%         cx = mod(cx,Lx);
%         cy = mod(cy,Ly);
%         quiver(cx,cy,lambda(1)*V(1,1),lambda(1)*V(2,1),'-w','linewidth',2);
%         quiver(cx,cy,lambda(2)*V(1,2),lambda(2)*V(2,2),'-w','linewidth',2);
    end
    
    % make colorbar
    if colorOpt == 1
        colormap(jet(NCLR));
        cb = colorbar;
        cb.Ticks = linspace(0,1,NCLRBAR+1);
        cb.TickLabels = cbTickCell;
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize = 18;
    end
        
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


%% Draw state with lowest coordination closest to z = 3

zmean = mean(zc,2);
dz = abs(zmean - 3);
[~, ff] = min(dz);

% reset figure for this frame
figure(20), clf, hold on, box on;

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

% plot box
plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
axis equal;
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.XLim = [-0.25 1.25]*L;
ax.YLim = [-0.25 1.25]*L;
