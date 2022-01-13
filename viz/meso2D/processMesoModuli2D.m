function processMesoModuli2D(floc,fpattern,savestr)
%% FUNCTION to analyze enthalpy-minimized mesophyll network simulations
% ASSUME ONLY HAS MODULI DATA FOR NOW, DESPITE .hess FILE FORMAT

% get list of files that fit pattern
flist = dir([floc '/' fpattern '*.posctc']);
NSIMS = length(flist);
if NSIMS == 0
    error('processMesoNetwork2D:noFilesFound','No files found with pattern %s in location %s\n',fpattern,floc);
else
    fprintf('Found %d file using pattern %s, processing...\n',NSIMS,fpattern);
end

% extract number of cells
NCELLS = sscanf(fpattern,'mesoDM2D_N%f_');

%% Loop over files, save information

% files to skip
fskip = false(NSIMS,1);

% data to save
fnameList = cell(NSIMS,1);
NFRAMESList = zeros(NSIMS,1);

nvList = cell(NSIMS,1);
LList = cell(NSIMS,1);
phiList = cell(NSIMS,1);
calAList = cell(NSIMS,1);
calA0List = cell(NSIMS,1);
SxxList = cell(NSIMS,1);
SyyList = cell(NSIMS,1);
SxyList = cell(NSIMS,1);
zList = cell(NSIMS,1);
t0List = cell(NSIMS,1);
kbList = cell(NSIMS,1);
cxList = cell(NSIMS,1);
cyList = cell(NSIMS,1);
BList = cell(NSIMS,1);
GList = cell(NSIMS,1);
prList = cell(NSIMS,1);

% Loop
for ss = 1:NSIMS
    % file info
    fname = flist(ss).name;
    fstr = [floc '/' fname];
    finfo = dir(fstr);
    if finfo.bytes == 0
        fprintf('** %s is empty, skipping...\n',fname);
        fskip(ss) = true;
        continue;
    end
    mesoData = readMesoNetworkCTCS2D(fstr);
    NFRAMES = mesoData.NFRAMES;
    NFRAMES = NFRAMES - 1;
    
    % also read in data from .hess (ASSUME THAT IT ONLY HAS MODULI DATA)
    hessname = [fname(1:end-7) '.hess'];
    hessstr = [floc '/' hessname];
    
    % Read in Hessian data
    fid = fopen(hessstr);

    % loop over frames, save moduli +  eigenvalues
    G = zeros(NFRAMES,1);
    B = zeros(NFRAMES,1);
%     mvals = cell(NFRAMES,1);
%     hvals = cell(NFRAMES,1);
%     svals = cell(NFRAMES,1);
    fprintf('** Reading in data from Moduli file %s\n',hessname);
    for ff = 1:NFRAMES
        % read in newfr
        newfrstr = fgetl(fid);

        % read in packing fraction
        phitmp = textscan(fid,'PACKF %f',1);

        % read in box size
        Ltmp = textscan(fid,'BOXSZ %f %f',1);
        emptystr = fgetl(fid);

        % read in shear modulus
        Gtmp = textscan(fid,'SHRMD %f',1);
        emptystr = fgetl(fid);
        G(ff) = Gtmp{1};

        % read in bulk modulus
        Btmp = textscan(fid,'BLKMD %f',1);
        emptystr = fgetl(fid);
        B(ff) = Btmp{1};

    %     % read in dynamical matrix eigenvalues
    %     mevalsstr = fgetl(fid);
    %     mvals{ff} = sscanf(mevalsstr(6:end),'%f');
    %     NMVALS = length(mvals{ff});
    %     fprintf('%s, %d mvals found\n',mevalsstr(1:5),NMVALS);
    %     
    %     hevalsstr = fgetl(fid);
    %     hvals{ff} = sscanf(hevalsstr(6:end),'%f');
    %     NHVALS = length(hvals{ff});
    %     fprintf('%s, %d hvals found\n',hevalsstr(1:5),NHVALS);
    %     
    %     sevalsstr = fgetl(fid);
    %     svals{ff} = sscanf(sevalsstr(6:end),'%f');
    %     NSVALS = length(svals{ff});
    %     fprintf('%s, %d hvals found\n',sevalsstr(1:5),NSVALS);

        % read in endfr
        endfrstr = fgetl(fid);
    end
    poissonRatio = (B-G)./(B+G);
    
    % save moduli info
    GList{ss} = G(2:NFRAMES);
    BList{ss} = B(2:NFRAMES);
    prList{ss} = poissonRatio(2:NFRAMES);
    
    % sim info
    NCELLS = mesoData.NCELLS;
    nv = mesoData.nv(2:NFRAMES,:);
    L = mesoData.L(2:NFRAMES,:);
    x = mesoData.x(2:NFRAMES,:);
    y = mesoData.y(2:NFRAMES,:);
    zc = mesoData.zc(2:NFRAMES,:);
    a0 = mesoData.a0(2:NFRAMES,:);
    l0 = mesoData.l0(2:NFRAMES,:);
    t0 = mesoData.t0(2:NFRAMES,:);
    kb = mesoData.kb(2:NFRAMES,:);
    S = mesoData.S(2:NFRAMES,:);
    phi = mesoData.phi(2:NFRAMES);
    p = mesoData.p(2:NFRAMES,:);
    a = mesoData.a(2:NFRAMES,:);
    
    NFRAMES = NFRAMES - 1;
    if NFRAMES == 0
        fprintf('** %s has 0 frames, skipping...\n',fname);
        fskip(ss) = true;
        continue;
    end
    NFRAMESList(ss) = NFRAMES;
    fnameList{ss} = fname;
    
    % cell vertices
    nvList{ss} = nv;
    
    % box
    LList{ss} = L;
    
    % preferred shape / cx
    calA0 = zeros(NFRAMES,NCELLS);
    cx = zeros(NFRAMES,NCELLS);
    cy = zeros(NFRAMES,NCELLS);
    for ff = 1:NFRAMES
        a0tmp = a0(ff,:);
        l0tmp = l0(ff,:);
        for cc = 1:NCELLS
            p0tmp = sum(l0tmp{cc});
            calA0(ff,cc) = p0tmp^2/(4.0*pi*a0tmp(cc));
        end
        xtmp = x(ff,:);
        ytmp = y(ff,:);
        cx(ff,:) = cellfun(@mean,xtmp);
        cy(ff,:) = cellfun(@mean,ytmp);
    end
    cxList{ss} = cx;
    cyList{ss} = cy;
    calA0List{ss} = calA0;

    % actual shape
    calAList{ss} = p.^2./(4.0*pi*a);
    
    % packing fraction
    phiList{ss} = phi;
    
    % connections
    zList{ss} = zc;
    
    % bending
    t0List{ss} = t0;
    kbList{ss} = kb;
    
    % stress
    SxxList{ss} = S(:,1);
    SyyList{ss} = S(:,2);
    SxyList{ss} = S(:,3);
    
    % sims to save
    svidx = ~fskip(1:ss);
    
    % save matfile as you go
    saveStruct.fnameList = fnameList(svidx);
    saveStruct.NFRAMESList = NFRAMESList(svidx);
    saveStruct.phiList = phiList(svidx,:);
    saveStruct.LList = LList(svidx,:);
    saveStruct.cxList = cxList(svidx,:);
    saveStruct.cyList = cyList(svidx,:);
    saveStruct.calA0List = calA0List(svidx,:);
    saveStruct.calAList = calAList(svidx,:);
    saveStruct.zList = zList(svidx,:);
    saveStruct.SxxList = SxxList(svidx);
    saveStruct.SxyList = SxyList(svidx);
    saveStruct.SyyList = SyyList(svidx);
    saveStruct.nvList = nvList(svidx,:);
    saveStruct.GList = GList(svidx);
    saveStruct.BList = BList(svidx);
    saveStruct.prList = prList(svidx);
    saveStruct.NSIMS = sum(svidx);
    saveStruct.NCELLS = NCELLS;
    save(savestr,'-struct','saveStruct');
end

% number of total, non-skipped sims
NSIMS = saveStruct.NSIMS;

%% Add ensemble averaged quantities to save struct

fprintf('Computing ensemble averages\n');

% number of frames
NFRAMESList = saveStruct.NFRAMESList;

% number of porosity bins
nporobins = 15;

% reload finalized data
phiList = saveStruct.phiList;
calAList = saveStruct.calAList;
LList = saveStruct.LList;
zList = saveStruct.zList;

% lists for later ensemble averaging
NFRAMETOT = sum(NFRAMESList);
poros = zeros(NFRAMETOT,1);
calAMeans = zeros(NFRAMETOT,1);
calAStds = zeros(NFRAMETOT,1);
zMeans = zeros(NFRAMETOT,1);
zStds = zeros(NFRAMETOT,1);
GVals = zeros(NFRAMETOT,1);
BVals = zeros(NFRAMETOT,1);
prVals = zeros(NFRAMETOT,1);


% construct binned porosity (skip first frame, always)
minPoro = 10;
maxPoro = 0;
last = 1;
for ss = 1:NSIMS
    % phi for sim ss
    phitmp = phiList{ss};
    
    % porosity
    porotmp = 1 - phitmp;
    
    % get min and max
    minporotmp = min(porotmp(1:end));
    maxporotmp = max(porotmp(1:end));
    if minporotmp < minPoro
        minPoro = minporotmp;
    end
    if maxporotmp > maxPoro
        maxPoro = maxporotmp;
    end
    
    % save quantities for ensemble averaging
    NFRAMEStmp = NFRAMESList(ss);
    next = last + NFRAMEStmp - 1;
    poros(last:next) = porotmp;
    
    % get cell shape information
    calAtmp = calAList{ss};
    calAMeans(last:next) = mean(calAtmp,2);
    calAStds(last:next) = std(calAtmp,0,2);
    
    % get contact information
    ztmp = zList{ss};
    zMeans(last:next) = mean(ztmp,2);
    zStds(last:next) = std(ztmp,0,2);
    
    % get moduli information
    Gtmp = GList{ss};
    GVals(last:next) = Gtmp;
    
    Btmp = BList{ss};
    BVals(last:next) = Btmp;
    
    prtmp = prList{ss};
    prVals(last:next) = prtmp;
    
    % update for next time
    last = next + 1;
end
poro_be = linspace(minPoro - 0.05*abs(minPoro),maxPoro + 0.05*abs(maxPoro),nporobins+1);
poro_bc = 0.5*(poro_be(2:end) + poro_be(1:end-1));

% ensemble averaged data
calA = zeros(nporobins,2);
z = zeros(nporobins,2);
G = zeros(nporobins,2);
B = zeros(nporobins,2);
pr = zeros(nporobins,2);

% loop over porosity bins, get all data
for bb = 1:nporobins
    % bin
    pbl = poro_be(bb);
    pbr = poro_be(bb+1);
    
    % get all data entries that fit in bin
    binidx = poros < pbr & poros > pbl;
    
    % save
    calA(bb,1) = mean(calAMeans(binidx));
    calA(bb,2) = mean(calAStds(binidx));
    
    z(bb,1) = mean(zMeans(binidx));
    z(bb,2) = mean(zStds(binidx));
    
    G(bb,1) = mean(GVals(binidx));
    G(bb,2) = std(GVals(binidx));
    
    B(bb,1) = mean(BVals(binidx));
    B(bb,2) = std(BVals(binidx));
    
    pr(bb,1) = mean(prVals(binidx));
    pr(bb,2) = std(prVals(binidx));
end

% save ensemble averages to save struct
saveStruct.poro_bc = poro_bc;
saveStruct.calA = calA;
saveStruct.z = z;
saveStruct.G = G;
saveStruct.B = B;
saveStruct.pr = pr;
save(savestr,'-struct','saveStruct');


fprintf('Wrote data to %s!\n',savestr);
fprintf('Ending.\n');


end





