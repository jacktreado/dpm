%% Draw dpm config files
% NOTE: need to get configs with contact info

clear;
close all;
clc;

% create file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb0.4_be100_h1_da0.5_dl5_cL0.5_cB4_t0m0.2_P1e-3_seed14.posctc';
% fstr = 'local/mesoDM2D_data/mesoDM2D_N32_n32_ca1.14_kl1_kb01e-3_be50_da0.02_dl10_P1e-4_seed27.posctc';
% fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.45;
NSKIP = 0;
idx(1:NSKIP) = zeros(NSKIP,1);
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
a = mesoData.a(idx,:);
a0 = mesoData.a0(idx,:);
l0 = mesoData.l0(idx,:);
t0 = mesoData.t0(idx,:);
kb = mesoData.kb(idx,:);
phi0 = sum(a0,2)./(LList(:,1).*LList(:,2));
phiA = sum(a,2)./(LList(:,1).*LList(:,2));

% loop over frames, get patch phi
phipatch = zeros(NFRAMES,1);
for ff = 1:NFRAMES
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    apatch = zeros(NCELLS,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        cx = mean(xtmp);
        cy = mean(ytmp);
        rx = xtmp - cx;
        ry = ytmp - cy;
        rads = sqrt(rx.^2 + ry.^2);
        xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
        ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
        apatch(nn) = polyarea(xtmp,ytmp);
    end
    phipatch(ff) = sum(apatch)/(LList(ff,1)*LList(ff,2));
end

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
P = mesoData.P(idx);
Sxx = mesoData.S(idx,1);
Syy = mesoData.S(idx,2);
Pvirial = 0.5*(Sxx + Syy);

% angle data
t0mean = cellfun(@mean,t0);
t0min = cellfun(@min,t0);
t0max = cellfun(@max,t0);

% print if multiple frames
if NFRAMES > 2
    
    figure(11), clf, hold on, box on;
    plot(1:NFRAMES,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(calA,2),'k-','linewidth',3);
    xlabel('frame id','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    figure(31), clf, hold on, box on;
    plot(1:NFRAMES,t0mean./pi,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,t0min./pi,':','color',[0.6 0.6 0.6],'linewidth',2);
    plot(1:NFRAMES,t0max./pi,'--','color',[0.25 0.25 0.25],'linewidth',1.5);
    xlabel('frame id','Interpreter','latex');
    ylabel('$\theta_{0k}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
   
    % plot contact network
    figure(21), clf, hold on, box on;
    plot(2:NFRAMES,zc(2:end,:),'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(2:NFRAMES,mean(zc(2:end,:),2),'k--','linewidth',2.5);
    xlabel('frame id.','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    % print pressure
    figure(12), clf, hold on, box on;
    plot(1:NFRAMES,P,'ko','markerfacecolor','b','markersize',10);
    plot(1:NFRAMES,Pvirial,'-k','linewidth',2);
    xlabel('frame id.','Interpreter','latex');
    ylabel('$P{\rm inst}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    legend({'$P = -(\partial U / \partial L)/(2L)$','$P_{\rm contacts}$'},'Interpreter','latex','fontsize',16,'location','best');
    
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
        plot(phipatch(2)-phipatch(2:end),zpop(2:end,cc),'-','linewidth',2,'color',clr(cc,:));
    end
    ylabel('$z$','Interpreter','latex');
    xlabel('$1-\phi$','Interpreter','latex');
    legend({'3','4','5','6','7','$\geq 8$'},'Interpreter','latex','fontsize',14,'location','best');
    ax = gca;
    ax.FontSize = 24;
    
    
    % plot area deviations
    figure(15), clf, hold on, box on;
    plot(1:NFRAMES,a0,'-','color',[0.5 0.5 0.5],'linewidth',1);
    plot(1:NFRAMES,mean(a0,2),'k-','linewidth',2.5);
    xlabel('frames','Interpreter','latex');
    ylabel('$a_0$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    
    % plot packing fraction and fraction perimeter at void
    figure(16), clf, hold on, box on;
    yyaxis right
    plot(1:NFRAMES-1,P(2:NFRAMES),'ks','markersize',10,'markerfacecolor','r');
    h = ylabel('$P$','Interpreter','latex');
    h.Color = 'r';
    ax = gca;
    ax.YColor = 'r';
    ax.YLim = [0 2*mean(P)];
    yyaxis left
    plot(1:NFRAMES-1,phipatch(2:NFRAMES),'ko','markersize',10,'markerfacecolor','b');
    h = ylabel('$\phi$','Interpreter','latex');
    h.Color = 'b';
    ax = gca;
    ax.YColor = 'b';
    ax.YLim = [0 1];
    xlabel('frame','Interpreter','latex');
    ax.FontSize = 24;
    ax.XLim = [1 NFRAMES];
    
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
    
    ct01CalAMean = totalData.ct01CalAMean;
    ct01CalAStd = totalData.ct01CalAStd;
    ct01Porosity = totalData.ct01Porosity;
    
    ct17CalAMean = totalData.ct17CalAMean;
    ct17CalAStd = totalData.ct17CalAStd;
    ct17Porosity = totalData.ct17Porosity;
    
    ct65CalAMean = totalData.ct65CalAMean;
    ct65CalAStd = totalData.ct65CalAStd;
    ct65Porosity = totalData.ct65Porosity;
    
    figure(17), clf, hold on, box on;
    errorbar(max(phipatch) - phipatch,mean(calA,2),std(calA,0,2),'-k','linewidth',1.5);
    errorbar(porosity,calAMean,calAMin,calAMax,'k>','markersize',10,'markerfacecolor','b');
    errorbar(porosity,calAMeas,dCalAMin,dCalAMax,'k>','markersize',14,'markerfacecolor','r');
    errorbar(ct01Porosity,ct01CalAMean,ct01CalAStd,'ko','markersize',12,'markerfacecolor','b');
    errorbar(ct17Porosity,ct17CalAMean,ct17CalAStd,'ks','markersize',12,'markerfacecolor','g');
    errorbar(ct65Porosity,ct65CalAMean,ct65CalAStd,'kd','markersize',12,'markerfacecolor','r');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    xlabel('$\phi_{\rm max}-\phi$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 24;
    ax.XLim = [0 0.2];
    ax.YLim = [1 1.5];
    
    porosity = ambroseData.porosity;
    aMean = ambroseData.areaMeanY;
    aSim = mean(a,2);
    a0Sim = mean(a0,2);
    
    figure(18), clf, hold on, box on;
    plot(porosity,aMean./aMean(1),'bo','markersize',10,'markerfacecolor','b');
    plot(max(phiA) - phiA,aSim./aSim(2),'ks','markersize',10);
    plot(max(phiA) - phiA,a0Sim./a0Sim(2),'k^','markersize',10);
    ylabel('$a/a(0)$','Interpreter','latex');
    xlabel('$1-\phi$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
    ax.YScale = 'log';
    legend({'data','$a/a(0)$','$a_0/a_0(0)$'},'Interpreter','latex','Fontsize',18);
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
colorOpt = 3;

if colorOpt == 1
    % color by real shape
    NCLR = 100;
%     calABins = linspace(0.999*min(calA(:)),1.001*max(calA(:)),NCLR+1);
    calABins = linspace(0.99,3,NCLR-1);
    calABins = [calABins 10000];
    cellCLR = jet(NCLR);
elseif colorOpt == 2
    % color by preferred shape
    NCLR = 100;
    calA0Bins = linspace(0.999*min(calA0(:)),1.001*max(calA0(:)),NCLR+1);
    cellCLR = jet(NCLR);
elseif colorOpt == 3
    % color face as gray
    t0_face_color = [0.7 0.7 0.7];
    
    % get bins for preferred curvature
    NCLRS = 100;
    t0_bins = linspace(-0.5*pi,0.5*pi,NCLRS-1);
    t0_bins = [-1e5 t0_bins 1e5];
    t0_clr_list = jet(NCLRS);
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
    FSTART = 21;
    FEND = FSTART;
%     FEND = NFRAMES;

    % set step size
    FSTEP = 1;
    DF = FEND - FSTART;
    if DF > 80 && DF <= 400
        FSTEP = 5;
    elseif DF > 400
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
    moviestr = 'mesoHMin2D_N32_n32_ca1.14_kb0.4_be200_h1_da0.5_dl5_cL0.5_cB4_t0m0.5_P1e-6_seed11.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
    ctccopy = -1:1;
end

% construct dl array for coloring
dl = cell(NFRAMES,NCELLS);
for ff = 1:NFRAMES
    for nn = 1:NCELLS
        dl{ff,nn} = zeros(nv(ff,nn),1);
    end
end

fnum = 1;
figure(fnum), clf, hold on, box on;
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    fprintf('printing frame ff = %d/%d, phi=%0.3g, phi0=%0.3g\n',ff,FEND,phipatch(ff),phi0(ff));
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    t0f = t0(ff,:);
    zctmp = zc(ff,:);
    zvff = zv(ff,:);
    L = LList(ff,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        t0tmp = t0f{nn};
        cx = mean(xtmp);
        cy = mean(ytmp);
        rtmp = rf{nn};
        nvtmp = nv(ff,nn);
        zvtmp = zvff{nn};
        
        % get color info
        switch colorOpt
            case 1
                cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
                clr = cellCLR(cbin,:);
            case 2
                cbin = calA0(ff,nn) > calA0Bins(1:end-1) & calA0(ff,nn) < calA0Bins(2:end);
                clr = cellCLR(cbin,:);
            case 3
                clr = t0_face_color;
            otherwise
                clr = cellCLR(IC(ff,nn),:);
        end
        if showverts == 1 && (nn == 29 || nn == 6)
            vpos = [xtmp, ytmp];
            patch('Faces',[1:nvtmp 1],'vertices',vpos,'FaceColor','none','EdgeColor','k','Linewidth',3,'LineStyle','-');
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if zvtmp(vv) > 0
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',2);
                        else
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',2);
                        end
                    end
                end
            end
        elseif showverts == 0
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
            ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
            for xx = -1:1
                for yy = -1:1
                    vpos = [xtmp + xx*L, ytmp + yy*L];
                    finfo = [1:nvtmp 1];
                    if colorOpt == 3
                        t0_clr = zeros(nvtmp,3);
                        for vv = 1:nvtmp
                            t0_bin = t0tmp(vv) > t0_bins(1:end-1) & t0tmp(vv) < t0_bins(2:end);
                            t0_clr(vv,:) = t0_clr_list(t0_bin,:);
                        end
                        patch('Faces',finfo,'vertices',vpos,'FaceVertexCData',t0_clr,'FaceColor',clr,'EdgeColor','interp','Linewidth',3);
                    else
                        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',1.5,'markersize',10);
                    end
%                     text(cx,cy,num2str(nn));
                end
            end
        end
    end
    
%     % plot vv contacts
%     gijtmp = gijList{ff};
%     xa = cell2mat(xf');
%     ya = cell2mat(yf');
%     NVTOT = sum(nv(ff,:));
%     for gi = 1:NVTOT
%         xi = xa(gi);
%         yi = ya(gi);
%         for gj = (gi+1):NVTOT
%             if (gijtmp(gi,gj) == 1)
%                 dx = xa(gj) - xi;
%                 dx = dx - L*round(dx/L);
%                 dy = ya(gj) - yi;
%                 dy = dy - L*round(dy/L);
%                 for xx = ctccopy
%                     for yy = ctccopy
%                         plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k-','linewidth',1.5);
%                     end
%                 end
%             end
%         end
%     end
        
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

return;
%% Draw individual cell shape dynamics

% pick random cell
cellidx = 5;

% reset figure for this frame
figure(20), clf, hold on, box on;

% get geometric info
xf = x(:,cellidx);
yf = y(:,cellidx);
nvf = nv(:,cellidx);
zvf = zv(:,cellidx);
% frList = [1 round(0.05*NFRAMES) round(0.1*NFRAMES) round(0.15*NFRAMES) round(0.2*NFRAMES) round(0.25*NFRAMES) round(0.6*NFRAMES) NFRAMES];
frList = [3:19];
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
    ytmp = ytmp - cy;
%     ytmp = ytmp - cy - (ii-1)*0.175*LP;
    ip1 = [2:nvtmp 1];
    for vv = 1:nvtmp
        if zvtmp(vv) <= 0
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',3,'color','r');
        else
            plot([xtmp(vv) xtmp(ip1(vv))],[ytmp(vv) ytmp(ip1(vv))],'-','linewidth',5,'color','k');
        end
    end
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = [-0.15*LP 0.15*LP];
    ax.YLim = [-0.15*LP 0.15*LP];
%     ax.YLim = [-NFR*0.175*LP 0.1*LP];
end
