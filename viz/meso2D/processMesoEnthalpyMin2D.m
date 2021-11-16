function processMesoEnthalpyMin2D(floc,fpattern,savestr)
%% FUNCTION to analyze enthalpy-minimized mesophyll network simulations
% NOTE: skip the first frame by default, not saved in ensemble or list 

% get list of files that fit pattern
flist = dir([floc '/' fpattern '*.posctc']);
NSIMS = length(flist);
if NSIMS == 0
    error('processMesoNetwork2D:noFilesFound','No files found with pattern %s in location %s\n',fpattern,floc);
else
    fprintf('Found %d file using pattern %s, processing...\n',NSIMS,fpattern);
end

% extract number of cells
NCELLS = sscanf(fpattern,'mesoHMin2D_N%f_');

%% Loop over files, save information

% files to skip
fskip = false(NSIMS,1);

% data to save
fnameList = cell(NSIMS,1);
NFRAMESList = zeros(NSIMS,1);

nvList = cell(NSIMS,NCELLS);
LList = cell(NSIMS,1);
phiList = cell(NSIMS,1);
polyList = cell(NSIMS,1);
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
    if finfo.bytes == 0
        fprintf('** %s is empty, skipping...\n',fname);
        fskip(ss) = true;
        continue;
    end
    mesoData = readMesoNetworkCTCS2D(fstr);
    
    % find void polygons, get number of OK frames
    polyListTmp = voidPolys(mesoData);
    NFRAMES = size(polyListTmp,1);

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
    polyList{ss} = polyListTmp(2:NFRAMES);
    
    NFRAMES = NFRAMES - 1;
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
    saveStruct.polyList = polyList(svidx,:);
    saveStruct.calA0List = calA0List(svidx,:);
    saveStruct.calAList = calAList(svidx,:);
    saveStruct.zList = zList(svidx,:);
    saveStruct.SxxList = SxxList(svidx,:);
    saveStruct.SxyList = SxyList(svidx,:);
    saveStruct.SyyList = SyyList(svidx,:);
    saveStruct.nvList = nvList(svidx,:);
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

% polygon type info
minpoly = 3;
maxpoly = 12;
polycheck = minpoly:maxpoly;
NPOLYCHECK = length(polycheck);

% reload finalized data
phiList = saveStruct.phiList;
calAList = saveStruct.calAList;
polyList = saveStruct.polyList;

% lists for later ensemble averaging
NFRAMETOT = sum(NFRAMESList);
poros = zeros(NFRAMETOT,1);
calAMeans = zeros(NFRAMETOT,1);
calAStds = zeros(NFRAMETOT,1);
npolys = zeros(NFRAMETOT,1);
polyTypeList = zeros(NFRAMETOT,NPOLYCHECK);
polyAreaList = zeros(NFRAMETOT,NPOLYCHECK);
zMeans = zeros(NFRAMETOT,1);
zStds = zeros(NFRAMETOT,1);


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
    minporotmp = min(porotmp(2:end));
    maxporotmp = max(porotmp(2:end));
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
    
    calAtmp = calAList{ss};
    calAMeans(last:next) = mean(calAtmp,2);
    calAStds(last:next) = std(calAtmp,0,2);
    
    % get polygon information
    LListtmp = LList{ss};
    polytmp = polyList{ss};
    npolys(last:next) = cellfun(@length,polytmp);
    polytypes = zeros(NFRAMEStmp,NPOLYCHECK);
    polyAreaPop = zeros(NFRAMEStmp,NPOLYCHECK);
    for ff = 1:NFRAMEStmp
        polyconfigtmp = polytmp{ff};
        NPOLYS = length(polyconfigtmp);
        polytypetmp = cellfun(@length,polyconfigtmp);
        for pp = minpoly:(maxpoly-1)
            polytypes(ff,pp-(minpoly-1)) = sum(polytypetmp == pp);
        end
        polytypes(ff,end) = sum(polytypetmp > maxpoly);
        Abox = LListtmp(ff,1)*LListtmp(ff,2);
        for ii = 1:NPOLYS
            ptmp = polyconfigtmp{ii};
            NE = size(ptmp,1);
            atmp = polyarea(ptmp(:,1),ptmp(:,2));
            if NE < maxpoly
                cc = NE-2;
            else
                cc = maxpoly-2;
            end
            polyAreaPop(ff,cc) = polyAreaPop(ff,cc) + (atmp/Abox);
        end
    end
    polyTypeList(last:next,:) = polytypes;
    polyAreaList(last:next,:) = polyAreaPop;
    
    % get contact information
    ztmp = zList{ss};
    zMeans(last:next) = mean(ztmp,2);
    zStds(last:next) = std(ztmp,0,2);
    
    % update for next time
    last = next + 1;
end
poro_be = linspace(minPoro - 0.05*abs(minPoro),maxPoro + 0.05*abs(maxPoro),nporobins+1);
poro_bc = 0.5*(poro_be(2:end) + poro_be(1:end-1));

% ensemble averaged data
calA = zeros(nporobins,2);
z = zeros(nporobins,2);
npoly = zeros(nporobins,2);
polytypeMean = zeros(nporobins,NPOLYCHECK);
polytypeStd = zeros(nporobins,NPOLYCHECK);
polyAreaMean = zeros(nporobins,NPOLYCHECK);
polyAreaStd = zeros(nporobins,NPOLYCHECK);

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
    
    npoly(bb,1) = mean(npolys(binidx));
    npoly(bb,2) = std(npolys(binidx));
    
    polytypeMean(bb,:) = mean(polyTypeList(binidx,:),1);
    polytypeStd(bb,:) = std(polyTypeList(binidx,:),0,1);
    
    polyAreaMean(bb,:) = mean(polyAreaList(binidx,:),1);
    polyAreaStd(bb,:) = std(polyAreaList(binidx,:),0,1);
end

% save ensemble averages to save struct
saveStruct.poro_bc = poro_bc;
saveStruct.calA = calA;
saveStruct.z = z;
saveStruct.npoly = npoly;
saveStruct.polytypeMean = polytypeMean;
saveStruct.polytypeStd = polytypeStd;
saveStruct.polyAreaMean = polyAreaMean;
saveStruct.polyAreaStd = polyAreaStd;
save(savestr,'-struct','saveStruct');


fprintf('Wrote data to %s!\n',savestr);
fprintf('Ending.\n');


end




%% FUNCTION to get void polygons based on data from simulation

function polys = voidPolys(mesoData)

% packing fraction (only take frames with phi > 0.25)
phi = mesoData.phi;
idx = phi > 0.01;

% number of frames
NFRAMES = sum(idx);

% sim info
NCELLS = mesoData.NCELLS;
nv = mesoData.nv(idx,:);
LList = mesoData.L;
ctcList = mesoData.ctcs;
x = mesoData.x(idx,:);
y = mesoData.y(idx,:);

% construct list of vv contacts
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
polys = cell(NFRAMES,1);
NPOLYS = zeros(NFRAMES,1);
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
    
    % get polygons
    cx = cellfun(@mean,x(ff,:))';
    cy = cellfun(@mean,y(ff,:))';
    L = LList(ff,1);
    [mainTiling, NFMAIN] = getMesoVoidPolygons(cx,cy,cijtmp,L);
    polys{ff} = mainTiling;
    NPOLYS(ff) = size(polys{ff},1);
    if NPOLYS(ff) >= 5
        fprintf('* On frame %d / %d, got %d void polygons...\n',ff,NFRAMES,NFMAIN);
    else
        fprintf('* On frame %d / %d, got %d void polygons, so ending here...\n',ff,NFRAMES,NFMAIN);
        polys(ff+1:end) = [];
        break;
    end
end

end
