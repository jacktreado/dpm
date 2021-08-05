function processCrawlingMonolayer(N, floc, saveloc, fpattern)
%% FUNCTION to process data from crawling monolayer simulations
% test: processCrawlingMonolayer(64, 'local/pos', '', 'mono_NT2e6_N64_n24_tc1.2_l10.005_Dr0.01_gamtt0_tau1000'); 
% Assuming fpattern defines an ensemble
% Assuming all sims have N particles

% get list of files in floc that fit pattern
fprintf('** In processCrawlingMonolayer\n');
fprintf('** Searching for files with pattern %s in %s \n',fpattern,floc);
if strcmp(floc(end),'/')
    floc(end) = [];
end
flist = dir([floc '/' fpattern '*.pos']);
NF = length(flist);
if NF == 0
    error('processCrawlingMonolayer:noFilesFound','No files found that fit pattern %s in %s. Ending here.\n',fpattern,floc);
else
    fprintf('** Found %d files in %s that fit pattern \n',NF,floc);
end

%% Loop over files in floc, extract data for ensemble averages

% files to skip
fskip   = false(NF,1);

% data to save
NFRAMES = zeros(NF,1);
dt      = zeros(NF,1);
L       = zeros(NF,2);
n       = zeros(NF,N);

simname = cell(NF,1);
t       = cell(NF,1);
calA    = cell(NF,N);
psi     = cell(NF,N);
cx      = cell(NF,N);
cy      = cell(NF,N);
zc      = cell(NF,N);
zv      = cell(NF,N);

% loop over files, extract data
for ff = 1:NF
    % get file
    fname = flist(ff).name;
    fstr = [floc '/' fname];
    finfo = dir(fstr);
    fsize = finfo.bytes;
    
    % check file size
    if fsize == 0
        fprintf('** File %s is empty, skipping...\n',fname);
        fskip(ff) = true;
        continue;
    end
    
    % save file name for parameters
    simname{ff} = fname;
    
    % load in data
    monoConfigData = readTumor2D(fstr);
    
    % get number of frames
    nftmp = monoConfigData.NFRAMES;
    NFRAMES(ff) = nftmp;

    % parse sim info
    dt(ff) = monoConfigData.dt(1);
    t{ff} = monoConfigData.t;
    n(ff,:) = monoConfigData.nv(1,:);
    L(ff,:) = monoConfigData.L(1,:);
    
    zctmp = monoConfigData.zc;
    zvtmp = monoConfigData.zv;
    atmp = monoConfigData.a;
    ptmp = monoConfigData.p;
    psitmp = monoConfigData.psi;
    xtmp = monoConfigData.x;
    ytmp = monoConfigData.y;
    
    % shape parameter
    calAtmp = ptmp.^2./(4.0*pi*atmp);
    
    % save particle-based sim info
    for nn = 1:N
        % per cell
        zc{ff,nn} = zctmp(:,nn);
        zv{ff,nn} = zvtmp(:,nn);
        calA{ff,nn} = calAtmp(:,nn);
        psi{ff,nn} = psitmp(:,nn);
        
        % cell c.o.m. over time
        xvtmp = xtmp(:,nn);
        yvtmp = ytmp(:,nn);
        cx{ff,nn} = cellfun(@mean,xvtmp);
        cy{ff,nn} = cellfun(@mean,yvtmp);
    end
end

% excise missing data
NFRAMES(fskip) = [];
dt(fskip) = [];
L(fskip,:) = [];
n(fskip,:) = [];

simname(fskip) = [];
t(fskip) = [];
calA(fskip,:) = [];
psi(fskip,:) = [];
cx(fskip,:) = [];
cy(fskip,:) = [];
zc(fskip,:) = [];
zv(fskip,:) = [];

NF = sum(~fskip);

%% Compute correlation functions

% time bins
tbins = 100;
tmax = max(cell2mat(t));
thalf = 0.5*tmax;
tbine = linspace(0.0,thalf,tbins+1);
tbinc = 0.5*(tbine(2:end) + tbine(1:end-1));

% correlations to save
msd = zeros(tbins,NF);
msd0 = cell(NF,1);
psicorr = zeros(tbins,NF);
zccorr = zeros(tbins,NF);
zvcorr = zeros(tbins,NF);
calAcorr = zeros(tbins,NF);

% loop over sims, compute correlations + bin
for ff = 1:NF
    % print to console
    fprintf('Computing 2-pt time correlations for %s\n',simname{ff});
    
    % number of frames
    nftmp = NFRAMES(ff);
    
    % get msd (average over cells, not time points)
    cxtmp = cell2mat(cx(ff,:));
    cytmp = cell2mat(cy(ff,:));
    
    % compute displacements for each cell
    msd0tmp = zeros(nftmp-1,1);
    msdtmp = zeros(nftmp-1,1);
    for nn = 1:N
        % get displacements
        dx = cxtmp(:,nn) - cxtmp(:,nn)';
        dy = cytmp(:,nn) - cytmp(:,nn)';
        dr = dx.^2 + dy.^2;
        
        % save square displacements from origin
        msd0tmp = msd0tmp + dr(1,2:end)';
        
        % compute entire msd
        for kk = 1:nftmp-1
            msdtmp(kk) = msdtmp(kk) + mean(diag(dr,kk));
        end
    end
    msd0tmp = msd0tmp./N;
    msdtmp = msdtmp./N;
    
    % save from origin
    msd0{ff} = msd0tmp;
    
    % bin true average
    ttmp = t{ff};
    dtframe = ttmp(2:end) - ttmp(1:end-1);
    dtframe = cumsum(dtframe);
    for bb = 1:tbins
        tbl = tbine(bb);
        tbr = tbine(bb+1);
        tinds = dtframe > tbl & dtframe < tbr;
        msd(bb,ff) = mean(msdtmp(tinds));
    end
    
    % temporary trajectorys
    psitmp = cell2mat(psi(ff,:));
    zctmp = cell2mat(zc(ff,:));
    zvtmp = cell2mat(zv(ff,:));
    calAtmp = cell2mat(calA(ff,:));
    
    % compute two-pt correlation functions using fft
    NT = length(ttmp);
    ntest = round(0.5*NT);
    psicorrtmp = zeros(ntest,1);
    zccorrtmp = zeros(ntest,1);
    zvcorrtmp = zeros(ntest,1);
    calAcorrtmp = zeros(ntest,1);
    for nn = 1:N
        % psi correlations
        f_psi = fft(psitmp(:,nn));
        fC_psi = (f_psi.*conj(f_psi))/NT;
        corr2 = ifft(fC_psi);
        corr2_rescale = (corr2 - mean(psitmp(:,nn))^2)./(mean(psitmp(:,nn).^2) - mean(psitmp(:,nn))^2);
        psicorrtmp = psicorrtmp + corr2_rescale(1:ntest)./N;
        
        % zc correlations
        f_zc = fft(zctmp(:,nn));
        fC_zc = (f_zc.*conj(f_zc))/NT;
        corr2 = ifft(fC_zc);
        corr2_rescale = (corr2 - mean(zctmp(:,nn))^2)./(mean(zctmp(:,nn).^2) - mean(zctmp(:,nn))^2);
        zccorrtmp = zccorrtmp + corr2_rescale(1:ntest)./N;
        
        % zv correlations
        f_zv = fft(zvtmp(:,nn));
        fC_zv = (f_zv.*conj(f_zv))/NT;
        corr2 = ifft(fC_zv);
        corr2_rescale = (corr2 - mean(zvtmp(:,nn))^2)./(mean(zvtmp(:,nn).^2) - mean(zvtmp(:,nn))^2);
        zvcorrtmp = zvcorrtmp + corr2_rescale(1:ntest)./N;
        
        % calA correlations
        f_calA = fft(calAtmp(:,nn));
        fC_calA = (f_calA.*conj(f_calA))/NT;
        corr2 = ifft(fC_calA);
        corr2_rescale = (corr2 - mean(calAtmp(:,nn))^2)./(mean(calAtmp(:,nn).^2) - mean(calAtmp(:,nn))^2);
        calAcorrtmp = calAcorrtmp + corr2_rescale(1:ntest)./N;
    end
    
    % save correlation functions
    dtcorr = dtframe(1:ntest);
    for bb = 1:tbins
        tbl = tbine(bb);
        tbr = tbine(bb+1);
        tinds = dtcorr > tbl & dtcorr < tbr;
        psicorr(bb,ff) = mean(psicorrtmp(tinds));
        zccorr(bb,ff) = mean(zccorrtmp(tinds));
        zvcorr(bb,ff) = mean(zvcorrtmp(tinds));
        calAcorr(bb,ff) = mean(calAcorrtmp(tinds));
    end
end




%% Save data

savestr = [saveloc '/' fpattern '.mat'];
save(savestr,'NFRAMES','dt','L','n','simname','t','calA','psi',...
    'cx','cy','zc','zv','NF','msd','msd0','psicorr','zccorr','zvcorr','calAcorr');




end