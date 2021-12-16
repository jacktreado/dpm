%% Draw dpm config files
% NOTE: need to get configs with contact info

clear;
close all;
clc;

% create file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb01e-3_be50_da0.05_dl7_P1e-4_h0.5_cL0_cB0_seed1.posctc';
% fstr = 'local/mesoDM2D_data/mesoDM2D_N64_n32_ca1.14_kb02e-2_be200_da0.05_dl5_P1e-4_h0.5_cL0_cB0_seed27.posctc';
% fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.35 & phi < 1.0;
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
if NFRAMES > 2
    Sxx = S(:,1);
    Syy = S(:,2);
    Sxy = S(:,3);
    
    figure(10), clf, hold on, box on;
    
    plot(find(Sxx<0),abs(Sxx(Sxx<0)),'ks','markersize',10,'MarkerFaceColor','r');
    plot(find(Sxx>0),Sxx(Sxx>0),'rs','markersize',10);
    
    plot(find(Syy<0),abs(Syy(Syy<0)),'ko','markersize',10,'MarkerFaceColor','b');
    plot(find(Syy>0),Syy(Syy>0),'bo','markersize',10);
    
    xlabel('frame id. ','Interpreter','latex');
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
    plot(1:NFRAMES,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(calA,2),'k-','linewidth',3);
    xlabel('frame id','Interpreter','latex');
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
    plot(1:NFRAMES,zc,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(zc,2),'k--','linewidth',2.5);
    xlabel('frame id.','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    % plot contact populations
    zchecklist = 3:8;
    NZCHECK = length(zchecklist);
    zpop = zeros(NFRAMES,NZCHECK);
    for ff = 1:NFRAMES
        zctmp = zc(ff,:);
        for cc = 1:NZCHECK
            zpop(ff,cc) = sum(zctmp == zchecklist(cc))/NCELLS;
        end
    end
    figure(50), clf, hold on, box on;
    clr = jet(NZCHECK);
    for cc = 1:NZCHECK
        plot(phi(2)-phi(2:end),zpop(2:end,cc),'-','linewidth',2,'color',clr(cc,:));
    end
    ylabel('$z$','Interpreter','latex');
    xlabel('$1-\phi$','Interpreter','latex');
    legend({'3','4','5','6','7','$\geq 8$'},'Interpreter','latex','fontsize',14,'location','best');
    ax = gca;
    ax.FontSize = 24;
    
    
    % plot area deviations
    ea = 0.5*(a./a0 - 1).^2;
    figure(15), clf, hold on, box on;
    plot(phi,ea,'-','color',[0.5 0.5 0.5],'linewidth',1);
    plot(phi,mean(ea,2),'k-','linewidth',2.5);
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$U_a$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    ax.YScale = 'log';
    
     % plot packing fractions
    figure(16), clf, hold on, box on;
    yyaxis left
    plot((2:NFRAMES) - 1,phi(2:NFRAMES),'ko','markersize',10,'markerfacecolor','b');
    h = ylabel('$\phi$','Interpreter','latex');
    h.Color = 'b';
    ax = gca;
    ax.YColor = 'b';
    yyaxis right
    plot((2:NFRAMES) - 1,P(2:NFRAMES),'ks','markersize',10,'markerfacecolor','r');
    h = ylabel('$P$','Interpreter','latex');
    h.Color = 'r';
    ax = gca;
    ax.YColor = 'r';
    xlabel('frame','Interpreter','latex');
    ax.YLim = [0 2e-4];
    ax.FontSize = 36;
    
    % plot shape vs porosity
    ambroseData = load('/Users/jacktreado/Jamming/Flowers/structure/plant/ambroseMesoCells/ambroseData.mat');
    totalData = load('/Users/jacktreado/Jamming/Flowers/structure/plant/ambroseMesoCells/totalData.mat');
    porosity = totalData.porosity;
    calAMean = totalData.calAMean;
    calAMin = totalData.calAMin;
    calAMax = totalData.calAMax;
    calAMeas = totalData.calAMeas;
    dCalAMin = totalData.dCalAMin;
    dCalAMax = totalData.dCalAMax;
    
    figure(17), clf, hold on, box on;
    errorbar(phi(2)-phi(2:end),mean(calA(2:end,:),2),std(calA(2:end,:),0,2),'-k','linewidth',1.75);
%     errorbar(porosity,calAMean,calAMin,calAMax,'k>','markersize',10,'markerfacecolor','b');
    errorbar(porosity,calAMeas,dCalAMin,dCalAMax,'k>','markersize',10,'markerfacecolor','r');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    xlabel('$1-\phi$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    porosity = ambroseData.porosity;
    aMean = ambroseData.areaMeanY;
    aSim = mean(a,2);
    a0Sim = mean(a0,2);
    
    figure(18), clf, hold on, box on;
    plot(porosity,aMean./aMean(1),'bo','markersize',10,'markerfacecolor','b');
    plot(phi(2)-phi(2:end),aSim(2:end)./aSim(2),'ks','markersize',10);
    plot(phi(2)-phi(2:end),a0Sim(2:end)./a0Sim(2),'k^','markersize',10);
    ylabel('$a/a(0)$','Interpreter','latex');
    xlabel('$1-\phi$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    ax.YScale = 'log';
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
colorOpt = 1;

if colorOpt == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
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

% construct list of contacts
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
end

% get frames to plot
if showverts == 0
    FSTART = 2;
    FSTEP = 1;
    if NFRAMES > 50
        FSTEP = 2;
    elseif NFRAMES > 150
        FSTEP = 10;
    end
%     FEND = NFRAMES;
    FEND = FSTART;
else
    FSTART = 1;
    FSTEP = 1;
    FEND = NFRAMES;
end

% make a movie
makeAMovie = 0;
ctccopy = 0;
if makeAMovie == 1
%     moviestr = 'debug.mp4';
    moviestr = 'mesoHMin2D_N64_n32_ca1.14_kb01e-3_be50_da0.05_dl7_P1e-4_h0.5_cL0_cB0_seed1.mp4';
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
            vpos = [xtmp, ytmp];
            patch('Faces',[1:nvtmp 1],'vertices',vpos,'FaceColor','none','EdgeColor','k','Linewidth',1.5,'LineStyle','-');
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if nn > 0
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',2);
                        else
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','w','LineWidth',2);
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
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',1.5,'markersize',10);
%                     text(cx,cy,num2str(nn));
                end
            end
        end
    end
    
    % plot vv contacts
    gijtmp = gijList{ff};
    xa = cell2mat(xf');
    ya = cell2mat(yf');
    NVTOT = sum(nv(ff,:));
    for gi = 1:NVTOT
        xi = xa(gi);
        yi = ya(gi);
        for gj = (gi+1):NVTOT
            if (gijtmp(gi,gj) == 1)
                dx = xa(gj) - xi;
                dx = dx - L*round(dx/L);
                dy = ya(gj) - yi;
                dy = dy - L*round(dy/L);
                for xx = ctccopy
                    for yy = ctccopy
                        plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k--','linewidth',2);
                    end
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

return;
%% Draw individual cell shape dynamics

% pick random cell
cellidx = 19;

% reset figure for this frame
figure(20), clf, hold on, box on;

% get geometric info
xf = x(:,cellidx);
yf = y(:,cellidx);
nvf = nv(:,cellidx);
zvf = zv(:,cellidx);
frList = [2 6 10 13];
NFR = length(frList);
LP = LList(end,1);
for ii = 1:NFR
    ff = frList(ii);
    xtmp = xf{ff};
    ytmp = yf{ff};
    zvtmp = zvf{ff};
    nvtmp = nvf(ff);
    
    cx = mean(xtmp);
    cy = mean(ytmp);
    xtmp = xtmp - cx;
    ytmp = ytmp - cy - (ii-1)*0.175*LP;
    ip1 = [2:nvtmp 1];
    for vv = 1:nvtmp
        if zvtmp(vv) <= 0
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',2,'color','g');
        else
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',2,'color','b');
        end
    end
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.15*LP 0.15*LP];
    ax.YLim = [-NFR*0.175*LP 0.1*LP];
end
