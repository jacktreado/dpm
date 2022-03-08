%% Draw dpm config files
% NOTE: need to get configs with contact info

clear;
close all;
clc;

% create file name
% fstr = 'local/mesoHMin2D_data/mesoHMin2D_N32_n32_ca1.14_kl1_kb01e-3_be50_da2e-2_dl5_P-1e-2_seed13.posctc';
% fstr = 'local/mesoDM2D_data/mesoDM2D_N32_n32_ca1.14_kl1_kb01e-3_be50_da0.02_dl10_P1e-4_seed27.posctc';
fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.45;
idx(1:3) = zeros(3,1);
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
% phiA = zeros(NFRAMES,1);
% for ff = 1:NFRAMES
%     phitmp = 0.0;
%     for nn = 1:NCELLS
%         phitmp = phitmp + a(ff,nn) + pi*mean(l0{ff,nn})^2*(0.5*nv(ff,nn)-1);
%     end
%     phitmp = phitmp/(LList(ff,1)*LList(ff,2));
%     phiA(ff) = phitmp;
% end


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

% print if multiple frames
if NFRAMES > 2
    
    figure(11), clf, hold on, box on;
    plot(1:NFRAMES,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    plot(1:NFRAMES,mean(calA,2),'k-','linewidth',3);
    xlabel('frame id','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
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
    
     % plot packing fraction and fraction perimeter at void
    figure(16), clf, hold on, box on;
    yyaxis left
    plot((2:NFRAMES) - 1,phi(2:NFRAMES),'ko','markersize',10,'markerfacecolor','b');
    h = ylabel('$\phi$','Interpreter','latex');
    h.Color = 'b';
    ax = gca;
    ax.YColor = 'b';
    ax.YLim = [0 1];
    yyaxis right
    plot((2:NFRAMES) - 1,phiA(2:NFRAMES),'ks','markersize',10,'markerfacecolor','r');
    h = ylabel('$\phi_{\rm area}$','Interpreter','latex');
    h.Color = 'r';
    ax = gca;
    ax.YColor = 'r';
    ax.YLim = [0 1];
    xlabel('frame','Interpreter','latex');
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
    
    ct01CalAMean = totalData.ct01CalAMean;
    ct01CalAStd = totalData.ct01CalAStd;
    ct01Porosity = totalData.ct01Porosity;
    
    ct17CalAMean = totalData.ct17CalAMean;
    ct17CalAStd = totalData.ct17CalAStd;
    ct17Porosity = totalData.ct17Porosity;
    
    figure(17), clf, hold on, box on;
    errorbar(max(phiA) - phiA,mean(calA,2),std(calA,0,2),'-k','linewidth',1.75);
%     errorbar(porosity,calAMean,calAMin,calAMax,'k>','markersize',10,'markerfacecolor','b');
    errorbar(porosity,calAMeas,dCalAMin,dCalAMax,'k>','markersize',10,'markerfacecolor','r');
    errorbar(ct01Porosity,ct01CalAMean,ct01CalAStd,'k>','markersize',10);
    errorbar(ct17Porosity,ct17CalAMean,ct17CalAStd,'k>','markersize',10);
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
    calABins = linspace(0.99,2.86,NCLR);
    calABins = [calABins 10];
    cellCLR = jet(NCLR);
elseif colorOpt == 2
    % color by preferred shape
    NCLR = 100;
    calA0Bins = linspace(0.999*min(calA0(:)),1.001*max(calA0(:)),NCLR+1);
    cellCLR = jet(NCLR);
elseif colorOpt == 3
    % color face as gray
    t0_face_color = [0.9 0.9 0.9];
    
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
    FSTART = 1;
    FSTEP = 1;
    if NFRAMES > 50
        FSTEP = 5;
    elseif NFRAMES > 150
        FSTEP = 20;
    end
    FEND = NFRAMES;
%     FEND = FSTART;
else
    FSTART = NFRAMES;
    FSTEP = 1;
    FEND = FSTART;
end

% make a movie
makeAMovie = 0;
ctccopy = 0;
if makeAMovie == 1
%     moviestr = 'debug.mp4';
    moviestr = 'negCurvature.mp4';
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
        if showverts == 1
            vpos = [xtmp, ytmp];
            patch('Faces',[1:nvtmp 1],'vertices',vpos,'FaceColor','none','EdgeColor','k','Linewidth',1.5,'LineStyle','-');
            for vv = 1:nvtmp
                rv = rtmp(vv);
                xplot = xtmp(vv) - rv;
                yplot = ytmp(vv) - rv;
                for xx = -1:1
                    for yy = -1:1
                        if zvtmp(vv) > 0
                            rectangle('Position',[xplot + xx*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor','k','LineWidth',2);
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
                        plot([xi, xi + dx] + xx*L,[yi, yi + dy] + yy*L,'k-','linewidth',1.5);
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
