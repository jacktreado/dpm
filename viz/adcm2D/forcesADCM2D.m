%% Script to draw cells and to check forces from force output

clear;
close all;
clc;

% position and force information
dpmloc = '/Users/jacktreado/Jamming/CellSim/dpm';
posfname = 'pos.test';
frcfname = 'measureFriction.frc';
posfstr = [dpmloc '/' posfname];
frcfstr = [dpmloc '/' frcfname];

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

if NCELLS < 10
    pbc = 0;
else
    pbc = 1;
end

% get force information
fid = fopen(frcfstr);
ff = 1;
frcDataList = cell(NFRAMES,1);
while ~feof(fid)
    % store new frame
    newfrstr = fgetl(fid);
    if ~strcmp(newfrstr,'NEWFR')
        error('Issue with counting frames, ending.');
    else
        frcdataframe = [];
    end
    
    % parse new line
    nxtfrstr = fgetl(fid);
    if strcmp(nxtfrstr,'ENDFR')
        ff = ff + 1;
        continue;
    else
        while ~strcmp(nxtfrstr,'ENDFR')
            % print
            fprintf(['On frame %d, nxtfrstr = ' nxtfrstr '\n'], ff);

            % parse "next frame string"
            frcdatatmp = sscanf(nxtfrstr,'%f %f %f %f %f %f %f');
            frcdatatmp = frcdatatmp';
            
            % store data
            frcdataframe = [frcdataframe; frcdatatmp];
            
            % get next
            nxtfrstr = fgetl(fid);
        end
        if strcmp(nxtfrstr,'ENDFR')
            frcDataList{ff} = frcdataframe;
            ff = ff + 1;
        end
    end
end


%% Draw cells with forces


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
% FSTART = 2135;
% if FSTART > NFRAMES
%     FSTART = NFRAMES;
% end
% FEND = FSTART;

% % movie frames
FSTART = 1;
FEND = NFRAMES;

% set step size
FSTEP = 1;
DF = FEND - FSTART;
if DF > 100 && DF <= 400
    FSTEP = 1;
elseif DF > 400 && DF <= 800
    FSTEP = 4;
elseif DF > 800
    FSTEP = 10;
end

% make a movie
makeAMovie = 0;
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
%             drawSSProjPoly(xtmp, ytmp, rtmp, clr, L, L, pbc);
            patch('Faces', [1:nvtmp, 1], 'Vertices', [xtmp ytmp], 'EdgeColor', 'k', 'FaceColor', clr);
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
    ax.XLim = [1.4585    1.9268];
    ax.YLim = [1.1241    1.5923];

        % plot forces
    if ~isempty(frcDataList{ff})
        % loop over forces, plot
        fdl = frcDataList{ff};
        NFRC = size(fdl,1);
        for gg = 1:NFRC
            fdltmp = fdl(gg,:);
            p = fdltmp([3 4]);
            del = fdltmp([5 6]);
            plot(p(1), p(2), 'ko','markersize', 6, 'markerfacecolor', 'b');
            plot(p(1) + del(1), p(2) + del(2), 'ko', 'markersize', 8, 'markerfacecolor', 'r');
            plot(p(1) + [0 del(1)], p(2) + [0 del(2)], 'k-', 'linewidth', 2);
        end
        test = 1;
    end
    
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