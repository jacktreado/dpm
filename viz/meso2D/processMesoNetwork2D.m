function processMesoNetwork2D(floc,fpattern,savestr)
%% FUNCTION to analyze mesophyll network simulations

% get list of files that fit pattern
flist = dir([floc '/' fpattern '*.pos']);
NSIMS = length(flist);
if NSIMS == 0
    error('processMesoNetwork2D:noFilesFound','No files found with pattern %s in location %s\n',fpattern,floc);
else
    fprintf('Found %d file using pattern %s, processing...\n',NSIMS,fpattern);
end

% extract number of cells
NCELLS = sscanf(fpattern,'meso2D_N%f_');

%% Loop over files, save information

% files to skip
fskip = false(NSIMS,1);

% data to save
fnameList = cell(NSIMS,1);
nvList = zeros(NSIMS,NCELLS);
LList = zeros(NSIMS,2);
NFRAMESList = zeros(NSIMS,1);

phiList = cell(NSIMS,1);
phi0List = cell(NSIMS,1);
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


% Loop
for ss = 1:NSIMS
    % file info
    fname = flist(ss).name;
    fstr = [floc '/' fname];
    finfo = dir(fstr);
    if finfo.size == 0
        fprintf('** %s is empty, skipping...\n',fname);
        fskip(ss) = true;
    end
    mesoData = readMesoNetwork2D(fstr);
    
    % get number of frames
    NFRAMES = mesoData.NFRAMES;
    NFRAMESList(ss) = NFRAMES;
    fnameList{ss} = fname;

    % sim info
    NCELLS = mesoData.NCELLS;
    nv = mesoData.nv(1,:);
    L = mesoData.L(1,:);
    Lx = L(1);
    Ly = L(2);
    x = mesoData.x;
    y = mesoData.y;
    zc = mesoData.zc;
    a0 = mesoData.a0;
    l0 = mesoData.l0;
    t0 = mesoData.t0;
    kb = mesoData.kb;
    
    % cell vertices
    nvList(ss,:) = nv;
    
    % box
    LList(ss,:) = L;
    
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
    p = mesoData.p;
    a = mesoData.a;
    calAList{ss} = p.^2./(4.0*pi*a);
    
    % packing fraction
    phiList{ss} = mesoData.phi;
    phi0List{ss} = sum(a0,2)/(Lx*Ly);
    
    % connections
    zList{ss} = zc;
    
    % bending
    t0List{ss} = t0;
    kbList{ss} = kb;
    
    % stress
    S = mesoData.S;
    SxxList{ss} = S(:,1);
    SyyList{ss} = S(:,2);
    SxyList{ss} = S(:,3);
    
    % save matfile as you go
    saveStruct.fnameList = fnameList(~fskip(1:ss));
    saveStruct.NFRAMESList = NFRAMESList(~fskip(1:ss));
    saveStruct.phiList = phiList(~fskip(1:ss),:);
    saveStruct.phi0List = phi0List(~fskip(1:ss),:);
    saveStruct.LList = LList(~fskip(1:ss),:);
    saveStruct.cxList = cxList(~fskip(1:ss),:);
    saveStruct.cyList = cyList(~fskip(1:ss),:);
    saveStruct.calA0List = calA0List(~fskip(1:ss),:);
    saveStruct.calAList = calAList(~fskip(1:ss),:);
    saveStruct.zList = zList(~fskip(1:ss),:);
    saveStruct.SxxList = SxxList(~fskip(1:ss),:);
    saveStruct.SxyList = SxyList(~fskip(1:ss),:);
    saveStruct.SyyList = SyyList(~fskip(1:ss),:);
    saveStruct.nvList = nvList(~fskip(1:ss),:);
    saveStruct.NSIMS = sum(~fskip(1:ss));
    saveStruct.NCELLS = NCELLS;
    save(savestr,'-struct','saveStruct');
end







end