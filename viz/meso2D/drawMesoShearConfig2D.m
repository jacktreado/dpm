%% Draw dpm config files
% NOTE: need to get configs with contact info

clear;
close all;
clc;

% create file name

% file name str
fstr = '~/Jamming/CellSim/dpm/pos.test';

% read in data
mesoData = readMesoShearConfig2D(fstr);

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.01;
phi = phi(idx);

% number of frames
NFRAMES = sum(idx);

% sim info
NCELLS = mesoData.NCELLS;
nv = mesoData.nv(idx,:);
LList = mesoData.L(idx,:);
gamma = mesoData.gamma(idx);
P = mesoData.P(idx,:);
S = mesoData.S(idx,:);
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

% print shear stress
dgamma = gamma(2) - gamma(1);

% get linear approx
m = (S(2)-S(1))/dgamma;
b = S(1) - m*gamma(1);
yfit = m*gamma + b;

figure(11), clf, hold on, box on;
plot(gamma,S,'-ks','markersize',10);
plot(gamma,yfit,'k-','linewidth',2);
xlabel('$\gamma$','Interpreter','latex');
ylabel('$\Sigma_{xy}$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;

G = -0.5*(S(3:end) - S(1:end-2))/dgamma;
figure(12), clf, hold on, box on;
plot(gamma(2:end-1),G,'-ks','markersize',10);
xlabel('$\gamma$','Interpreter','latex');
ylabel('$G$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;

%% Draw cells

% information for ellipses
th = 0.0:0.01:(2.0*pi);
NTHE = length(th);
ex = cos(th);
ey = sin(th);


% show vertices or not
showverts = 0;

% color by shape or size
colorOpt = 0;

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
    moviestr = 'meso2D_shear.mp4';
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
    L = LList(ff,1);
    gamtmp = gamma(ff);
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
                        rectangle('Position',[xplot + xx*L + gamtmp*yy*L, yplot + yy*L, 2.0*rv, 2.0*rv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr,'LineWidth',0.2);
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
                    vpos = [xtmp + xx*L + gamtmp*yy*L, ytmp + yy*L];
                    finfo = [1:nvtmp 1];
                    patch('Faces',finfo,'vertices',vpos,'FaceColor',clr,'EdgeColor','k');
                end
            end
        end
    end
        
    % plot box
    plot([0 1 1+gamtmp gamtmp 0]*L,[0 0 L L 0],'-k','linewidth',2);
    plot([0 1 1 0 0]*L,[0 0 L L 0],'--r','linewidth',1.2);
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    
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
