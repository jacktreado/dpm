%% Script to plot mesophyll triptych with PBCs, free boundary growth, and to trace out phi & pressure

clear;
close all;
clc;

% file name
fstr = 'local/mesoHMin2D_data/mesoHMin2D_N64_n32_ca1.14_kb0.4_be100_h1_da0.5_dl5_cL0.5_cB4_t0m0.2_P1e-3_seed14.posctc';

% input pressure (UPDATE BY HAND EACH TIME)
P0 = 1e-3;

% read in data
mesoData = readMesoNetworkCTCS2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.35;
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
P = mesoData.P(idx,:);

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

% particle shape data
p = mesoData.p(idx,:);
a = mesoData.a(idx,:);
calA = p.^2./(4.0*pi*a);


%% Draw triptych

% color by real shape
NCLR = 100;
calABins = linspace(0.99,3,NCLR-1);
calABins = [calABins 10000];
cellCLR = jet(NCLR);

% colorbar
NCB = 6;
cbtickcell = cell(NCB,1);
for cc = 1:NCB
    binloc = round(((cc-1)/NCB)*NCLR)+1;
    cbtickcell{cc} = sprintf('%0.3g',calABins(binloc+1));
end

% movie frames
FSTART = 2;
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

% make a movie
makeAMovie = 1;
ctccopy = 0;
if makeAMovie == 1
    moviestr = 'mesoHMin2D_N64_n32_ca1.14_kb0.4_be100_h1_da0.5_dl5_cL0.5_cB4_t0m0.2_P1e-3_seed14_trip.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 2;
    open(vobj);
    ctccopy = -1:1;
end

for ff = FSTART:FSTEP:FEND
    % print info
    fprintf('printing frame ff = %d/%d, phi=%0.3g\n',ff,FEND,phipatch(ff));
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    
    % box info
    L = LList(ff,1);
    
    % need to clear last frame
    figobj = figure(1); clf, hold on, box on;
    figobj.Color = 'w';
    
    % loop over cells
    for nn = 1:NCELLS
        % cell info
        xtmp = xf{nn};
        ytmp = yf{nn};
        cx = mean(xtmp);
        cy = mean(ytmp);
        rtmp = rf{nn};
        nvtmp = nv(ff,nn);
        
        % color
        cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
        clr = cellCLR(cbin,:);
        
        % subplot 1
        subplot(1,3,1), hold on, box on;
        
        rx = xtmp - cx;
        ry = ytmp - cy;
        rads = sqrt(rx.^2 + ry.^2);
        xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
        ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
        for xx = -1:1
            for yy = -1:1
                vpos = [xtmp + xx*L, ytmp + yy*L];
                finfo = [1:nvtmp 1];
                patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',1.5,'markersize',10);
            end
        end
        
        % same to subplot 2
        subplot(1,3,2), hold on, box on;
        
        rx = xtmp - cx;
        ry = ytmp - cy;
        rads = sqrt(rx.^2 + ry.^2);
        xtmp = xtmp + sin(pi/3)*rtmp.*(rx./rads);
        ytmp = ytmp + sin(pi/3)*rtmp.*(ry./rads);
        rx = xtmp - cx;
        ry = ytmp - cy;
        cx = mod(cx,L);
        cy = mod(cy,L);
        xtmp = rx + cx;
        ytmp = ry + cy;
        vpos = [xtmp, ytmp];
        finfo = [1:nvtmp 1];
        patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k','Linewidth',1.5,'markersize',10);
    end
    
    % scale subplot 1
    ax1 = subplot(1,3,1); hold on, box on;

    plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
    ax1.XTick = [];
    ax1.YTick = [];
    ax1.XLim = [-0.25 1.25]*L;
    ax1.YLim = [-0.25 1.25]*L;
    xlabel('Tissue frame','Interpreter','latex','Fontsize',24);
    axis square;

    
    colormap(jet);
    cb = colorbar(ax1,'Location','westoutside');
    cb.Ticks = linspace(0,1,NCB);
    cb.TickLabels = cbtickcell;
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 24;

    
    % scale subplot 2
    ax2 = subplot(1,3,2); hold on, box on;
    
    plot([0 L L 0 0], [0 0 L L 0], 'k-', 'linewidth', 1.5);
    ax2.XTick = [];
    ax2.YTick = [];
    ax2.XLim = [-0.25 1.25]*LList(end,1);
    ax2.YLim = [-0.25 1.25]*LList(end,2);
    axis square;
    xlabel('Lab frame','Interpreter','latex','Fontsize',24);
    
    % plot subplot 3
    ax3 = subplot(1,3,3); hold on, box on;
    
    yyaxis right
    plot(1:ff-1,P(2:ff)./P0,'ks','markersize',10,'markerfacecolor','r');
    h = ylabel('$P/P_0$','Interpreter','latex');
    h.Color = 'r';
    ax3.YColor = 'r';
    ax3.YLim = [0.9 1.1];
    
    yyaxis left
    plot(1:ff-1,phipatch(2:ff),'ko','markersize',10,'markerfacecolor','b');
    h = ylabel('$\phi$','Interpreter','latex');
    h.Color = 'b';
    ax3_r = gca;
    ax3_r.YColor = 'b';
    ax3_r.YLim = [0 1];
    xlabel('step \#','Interpreter','latex');
    ax3_r.FontSize = 24;
    ax3.XLim = [1 NFRAMES-1];
    axis square;
    
    drawnow;
    
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
