%% Script to draw outcome of adcm2D simulation

clear;
close all;
clc;

% position information
dpmloc = '/Users/jacktreado/Jamming/CellSim/dpm';
posfname = 'adcm2D_test.pos';
posfstr = [dpmloc '/' posfname];

% load into simulation data struct
simdata = readADCM2DTrajectory(posfstr);

% parse data
NFRAMES = simdata.NFRAMES;
NCELLS = simdata.NCELLS;
nv = simdata.nv;
LList = simdata.L;
S = simdata.S;
x = simdata.x;
y = simdata.y;
r = simdata.r;
p = simdata.p;
a = simdata.a;
st = simdata.st;
calA = p.^2 ./ (4.0 * pi .* a);
a0 = simdata.a0;
phi0 = sum(a0,2)./(LList(:,1).*LList(:,2));
phiA = sum(a,2)./(LList(:,1).*LList(:,2));

if NCELLS < 3
    pbc = 0;
else
    pbc = 1;
end

% stresses
P = S(:,1);
Sxy = S(:,2);

% phi = area + edge of circulolines
ravg = cellfun(@mean, r);
phiX = (sum(a, 2) + sum(ravg .* p, 2) + sum(pi .* ravg.^2, 2)) ./ (LList(:,1) .* LList(:,2));
figure(10), clf, hold on, box on;
plot(phiX, 'ko','markersize',8);
xlabel('frame','Interpreter','latex');
ylabel('$\phi$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;


% plot surface tension fluctuations
plotClr = winter(3);
meanst = cellfun(@mean, st);
figure(11), clf, hold on, box on;
patchErrorBar((1:size(meanst,1))',mean(meanst,2),std(meanst,0,2),'-o',8,'b','b',0.5);
xlabel('frame','Interpreter','latex');
ylabel('$\langle \gamma_{i\mu} \rangle$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;


% plot perimeter fluctuations
figure(12), clf, hold on, box on;
patchErrorBar((1:size(meanst,1))',mean(p - mean(p),2),std(p - mean(p),0,2),'-o',8,plotClr(1,:),plotClr(1,:),0.5);
% plot(p - mean(p), '-o','markersize',8);
xlabel('frame','Interpreter','latex');
ylabel('$\Delta p_\mu$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;

% plot area fluctuations
figure(13), clf, hold on, box on;
patchErrorBar((1:size(meanst,1))',mean(a - mean(a),2),std(a - mean(a),0,2),'-o',8,plotClr(1,:),plotClr(1,:),0.5);
% plot(p - mean(p), '-o','markersize',8);
xlabel('frame','Interpreter','latex');
ylabel('$\Delta a_\mu$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;


% plot area strains
da = a./a0 - 1.0;
figure(14), clf, hold on, box on;
patchErrorBar((1:NFRAMES)',mean(da,2),std(da,0,2),'-o',8,plotClr(1,:),plotClr(1,:),0.5);
% plot(p - mean(p), '-o','markersize',8);
xlabel('frame','Interpreter','latex');
ylabel('$\epsilon_a$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;


% plot shape fluctuations
figure(15), clf, hold on, box on;
patchErrorBar((1:size(meanst,1))',mean(calA,2),std(calA,0,2),'-o',8,plotClr(2,:),plotClr(2,:),0.5);
xlabel('frame','Interpreter','latex');
ylabel('$\langle \mathcal{A} \rangle$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;


% plot instantaneous pressure 
figure(16), clf, hold on, box on;
plot(1:NFRAMES,P,'-ko','markerfacecolor',plotClr(2,:),'markersize',10);
xlabel('frame','Interpreter','latex');
ylabel('$P$','Interpreter','latex');
ax = gca;
ax.FontSize = 24;
ax.YScale = 'log';


%% Draw cells over time

% color by shape or size
colorOpt = -1;

% color option
if colorOpt == 1
    % color by real shape
    NCLR = 100;
    calABins = linspace(0.99,3,NCLR-1);
    calABins = [calABins 10000];
    cellCLR = jet(NCLR);
elseif colorOpt == 2
    % color by area strain
    da = (a ./ a0) - 1.0;
    absda = abs(da);
    NCLR = 100;
    daBins = linspace(-1.01*max(absda(:)),1.01*max(absda(:)),NCLR+1);
    cellCLR = jet(NCLR);
else
    [nvUQ, ~, IC] = unique(nv);
    IC = reshape(IC,NFRAMES,NCELLS);
    NUQ = length(nvUQ);
    cellCLR = winter(NUQ);
end


% get frames to plot

% single frame
% FSTART = NFRAMES;
% if FSTART > NFRAMES
%     FSTART = NFRAMES;
% end
% FEND = FSTART;

% movie frames
FSTART = 1;
FEND = NFRAMES;

% set step size
FSTEP = 1;
DF = FEND - FSTART;
if DF > 100 && DF <= 400
    FSTEP = 2;
elseif DF > 400 && DF <= 800
    FSTEP = 6;
elseif DF > 800
    FSTEP = 10;
end

% make a movie
makeAMovie = 1;
if FSTART == FEND
    % draw boundary if single frame, don't make a movie
    bndry = 1;
    makeAMovie = 0;
else
    bndry = 0;
end
ctccopy = 1;
if makeAMovie == 1
    moviestr = 'adcm2D_VARVERTS_gam0.1_W0_l21.0_l10.01_dT1.5.mp4';
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
    fprintf('printing frame ff = %d/%d\n',ff,FEND);
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    L = LList(ff,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        rtmp = rf{nn};
        nvtmp = nv(ff,nn);
        
        % get color info
        switch colorOpt
            case -1
                clr = [0 0 1];
            case 1
                cbin = calA(ff,nn) > calABins(1:end-1) & calA(ff,nn) < calABins(2:end);
                clr = cellCLR(cbin,:);
            case 2
                cbin = da(ff,nn) > daBins(1:end-1) & da(ff,nn) < daBins(2:end);
                clr = cellCLR(cbin,:);
            otherwise
                clr = cellCLR(IC(ff,nn),:);
        end
        
        % draw SS polygon
        if bndry == 1
            drawSSPoly(xtmp, ytmp, rtmp, clr, L, L, pbc, bndry);
        else
            drawSSProjPoly(xtmp, ytmp, rtmp, clr, L, L, pbc);
        end
        text(mean(xtmp), mean(ytmp), num2str(nn));
    end
        
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


