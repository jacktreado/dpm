%% Draw dpm config files

clear;
close all;
clc;

% universal params
Nstr = '62';
nstr = '16';
calA0str = '1.14';
kbstr = '0';
seedstr = '10';

% simtype
simstr = 'a2j';

% info to edit
datadir = 'local/jam_data';

% anneal2Jam sims
if strcmp(simstr,'a2j')
    % a2j params
    trunstr = '50';
    T0str = '1e-2';
    
    % a2j file pattern
    fpattern = [simstr '_N' Nstr '_n' nstr '_calA0' calA0str '_kb' kbstr '_trun' trunstr '_T0' T0str '_seed' seedstr];
    
% lobed particle sims
elseif strcmp(simstr,'lobes')
    % a2j params
    thAstr = '1.0';
    thKstr = '1.0';
    
    % a2j file pattern
    fpattern = [simstr '_N' Nstr '_n' nstr '_calA0' calA0str 'kb' kbstr '_thA' thAstr '_thK' thKstr '_seed' seedstr];
else
    error('drawDPMConfig:simNotSupported','Error: %s simtype not yet supported, ending here. Check spelling and redo.\n',simstr);
end

% file info
vizdir = pwd;
fname = [fpattern '.pos'];
% fstr = [vizdir '/' datadir '/' fname];
fstr = '../../monodisperse.test';

% read in data
dpmData = readDPMConfig(fstr);

% get number of frames
NFRAMES = dpmData.NFRAMES;
NCELLS = dpmData.NCELLS;
nv = dpmData.nv(1,:);
LList = dpmData.L;
x = dpmData.x;
y = dpmData.y;
r = dpmData.r;
zc = dpmData.zc;
zv = dpmData.zv;
a0 = dpmData.a0;
l0 = dpmData.l0;

p = dpmData.p;
a = dpmData.a;
calA = p.^2./(4.0*pi*a);

S = dpmData.S;
P = 0.5*(S(:,1) + S(:,2));

phi = dpmData.phi;

%% Draw cells

% show vertices or not
showverts = 0;

% get cell colors
[nvUQ, ~, IC] = unique(nv);
NUQ = length(nvUQ);
allClr = jet(6);
cellCLR = allClr(1:NUQ,:);

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
    moviestr = 'compress2jamming_dpb.mp4';
    vobj = VideoWriter(moviestr,'MPEG-4');
    vobj.FrameRate = 15;
    open(vobj);
end

fnum = 1;
figure(fnum), clf, hold on, box on;
pactual = zeros(NFRAMES,NCELLS);
for ff = FSTART:FSTEP:FEND
    % reset figure for this frame
    figure(fnum), clf, hold on, box on;
    fprintf('printing frame ff = %d/%d\n',ff,FEND);
    
    % get geometric info
    xf = x(ff,:);
    yf = y(ff,:);
    rf = r(ff,:);
    zctmp = zc(ff,:);
    L = LList(ff,1);
    for nn = 1:NCELLS
        xtmp = xf{nn};
        ytmp = yf{nn};
        lx = xtmp([2:end 1]) - xtmp;
        ly = ytmp([2:end 1]) - ytmp;
        l = sqrt(lx.^2 + ly.^2);
        ptmp = sum(l);
        pactual(ff,nn) = ptmp;
        rtmp = rf{nn};
        clr = cellCLR(IC(nn),:);
%         clr = cellCLR(nn,:);
        if showverts == 1
            for vv = 1:nv(nn)
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
            cx = mean(xtmp);
            cy = mean(ytmp);
            rx = xtmp - cx;
            ry = ytmp - cy;
            rads = sqrt(rx.^2 + ry.^2);
            xtmp = xtmp + 0.8*rtmp.*(rx./rads);
            ytmp = ytmp + 0.8*rtmp.*(ry./rads);
            for xx = -1:1
                for yy = -1:1
                    vpos = [xtmp + xx*L, ytmp + yy*L];
                    finfo = [1:nv(nn) 1];
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
    
    plot(phi(Sxx<0),abs(Sxx(Sxx<0)),'ks','markersize',10,'MarkerFaceColor','r');
    plot(phi(Sxx>0),Sxx(Sxx>0),'ro','markersize',10);
    
    plot(phi(Syy<0),abs(Syy(Syy<0)),'ks','markersize',10,'MarkerFaceColor','b');
    plot(phi(Syy>0),Syy(Syy>0),'bo','markersize',10);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xx}$, $\Sigma_{yy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(11), clf, hold on, box on;
    
    plot(phi(Sxy<0),abs(Sxy(Sxy<0)),'ks','markersize',10,'MarkerFaceColor','k');
    plot(phi(Sxy>0),Sxy(Sxy>0),'ko','markersize',10);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\Sigma_{xy}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;
%     ax.YScale = 'log';
    
    figure(12), clf, hold on, box on;
    plot(phi,calA,'-','color',[0.5 0.5 0.5],'linewidth',1.2);
    errorbar(phi,mean(calA,2),std(calA,0,2),'k--','linewidth',2);
    
    xlabel('$\phi$','Interpreter','latex');
    ylabel('$\mathcal{A}$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 22;    
end