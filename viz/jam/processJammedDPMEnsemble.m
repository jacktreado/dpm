function processJammedDPMEnsemble(ensembleStr, savestr)
%% Process ensemble of jammed dpm configurations
% ** ASSUMING ONLY 1 FRAME, will skip configs that are empty or have
% multiple frames
% 08/11 To-do: use to process packing in individual svoro cells

% get files that make up ensemble
ensList = dir([ensembleStr '*.pos']);
NEN = length(ensList);
if NEN == 0
    error('processJammedDPMEnsemble:noFilesFound','\nNo files found using ensembleStr:%s. Ending here.\n',ensembleStr);
else
    fprintf('** Processing %d files stored in %s\n',NEN,ensembleStr);
end

%% Loop over files in ensemble, extract data

% files to skip
fskip       = false(NEN,1);

% data to save
filename            = cell(NEN,1);
NFRAMES             = zeros(NEN,1);
NCELLS              = zeros(NEN,1);        % number of cells in each sim
nv                  = cell(NEN,1);         % # of vertices on each particle
L                   = zeros(NEN,2);        % box lengths
phi                 = zeros(NEN,1);
S                   = zeros(NEN,3);
calA                = cell(NEN,1);
meanCalA            = zeros(NEN,1);
stdCalA             = zeros(NEN,1);
calA0               = cell(NEN,1);
abump               = cell(NEN,1);
zc                  = cell(NEN,1);
zv                  = cell(NEN,1);

% voronoi arrays
voroAreas           = cell(NEN,1);          % voronoi areas
voroCalA            = cell(NEN,1);          % voronoi calA
voroPackU           = cell(NEN,1);          % voronoi packing energy

% loop over ensemble
for ee = 1:NEN
    % get file
    fname = ensList(ee).name;
    floc = ensList(ee).folder;
    fstr = [floc '/' fname];
    finfo = dir(fstr);
    fsize = finfo.bytes;
    
    % check file size
    if fsize == 0
        fprintf('** File %s is empty, skipping...\n',fname);
        fskip(ee) = true;
        continue;
    end
    
    % save file name for parameters
    filename{ee} = fname;
    
    % load in data
    dpmConfigData = readDPMConfig(fstr);
    
    % check # of frames
    NFRAMEStmp = dpmConfigData.NFRAMES;
    if NFRAMEStmp ~= 1
        fprintf('** File %s has %d frames, skipping...\n',fname,NFRAMEStmp);
        fskip(ee) = true;
        continue;
    end
    NFRAMES(ee) = NFRAMEStmp;
    NCELLS(ee) = dpmConfigData.NCELLS;
    nv{ee} = dpmConfigData.nv;
    L(ee,:) = dpmConfigData.L;
    
    % save data
    phi(ee) = dpmConfigData.phi;
    
    Stmp    = dpmConfigData.S;
    S(ee,1) = Stmp(1);
    S(ee,2) = Stmp(2);
    S(ee,3) = Stmp(3);
    
    % shape info
    ptmp = dpmConfigData.p;
    atmp = dpmConfigData.a;
    calA{ee} = ptmp.^2./(4.0*pi*atmp);
    meanCalA(ee) = mean(calA{ee});
    stdCalA(ee) = std(calA{ee});
    
    % compute bumpy, exposed area
    r = dpmConfigData.r;
    abumptmp = zeros(NCELLS(ee),1);
    for nn = 1:NCELLS(ee)
        rn = r{nn}(1);
        nvn = nv{ee}(1);
        abumptmp(nn) = atmp(nn) + pi*(rn^2)*(0.5*nvn - 1);
    end
    abump{ee} = abumptmp;
    
    % compute calA0s
    l0tmp = dpmConfigData.l0;
    a0tmp = dpmConfigData.a0;
    calA0tmp = zeros(NCELLS(ee),1);
    for nn = 1:NCELLS(ee)
        p0tmp = sum(l0tmp{nn});
        calA0tmp(nn) = p0tmp^2/(4.0*pi*a0tmp(nn));
    end
    calA0{ee} = calA0tmp;
    
    % contact info
    zc{ee} = dpmConfigData.zc;
    zv{ee} = dpmConfigData.zv;
    
    % voronoi
    fprintf('* Computing Voronoi data for sim ...\n');
    x = dpmConfigData.x;
    y = dpmConfigData.y;
    [svoroFaceList, Vorig, voroCalATmp, voroAreasTmp] = getSurfaceVoronoi(x,y,nv{ee},L(ee,1));
    fprintf('...Voronoi done.\n');
    
    % save geometric info
    voroAreas{ee} = voroAreasTmp;
    voroCalA{ee} = voroCalATmp;
    
    % compute energies of all relaxed particles in voronoi cells
    Utmp = zeros(NCELLS(ee),1);
%     for nn = 1:NCELLS(ee)
%         % get info for voronoi cell around particle nn
%         finfo = svoroFaceList{nn};
%         V = Vorig(finfo,:);
%         
%         % pass to function to get energy
%         fprintf('* Getting relaxed single particle in cell %d...',nn);
%         Utmp(nn) = singleCellVoronoiPackingEnergy(nv{ee}(nn),x{nn},y{nn},r{nn},a0tmp(nn),l0tmp{nn}(1),V);
%         fprintf(' done.\n');
%     end
    voroPackU{ee} = Utmp;
    
    % save and append
    inds = ~fskip(1:ee);
    s.filename = filename(inds);
    s.NFRAMES = NFRAMES(inds);
    s.NFRAMES = NCELLS(inds);
    s.nv = nv(inds);
    s.L = L(inds,:);
    s.phi = phi(inds);
    s.S = S(inds);
    s.calA = calA(inds);
    s.meanCalA = meanCalA(inds);
    s.stdCalA = stdCalA(inds);
    s.calA0 = calA0(inds);
    s.abump = abump(inds);
    s.zv = zv(inds);
    s.zc = zc(inds);
    s.voroAreas = voroAreas(inds);
    s.voroCalA = voroCalA(inds);
    s.voroPackU = voroPackU(inds);
    sdata = whos('s');
    
    fprintf('On ee = %d, saving save file for struct of size %0.5g MB...\n',ee,sdata.bytes/1e6);
    tic;
    save(savestr,'-struct','s');
    toc;
end

% % delete extra entries
% fprintf('Delete unneeded entries ...\n');
% filename(fskip) = [];
% NFRAMES(fskip) = [];
% NCELLS(fskip) = [];
% nv(fskip) = [];
% L(fskip,:) = [];
% phi(fskip) = [];
% S(fskip,:) = [];
% calA(fskip) = [];
% meanCalA(fskip) = [];
% stdCalA(fskip) = [];
% calA0(fskip) = [];
% zv(fskip) = [];
% zc(fskip) = [];
% voroAreas(fskip) = [];
% voroCalA(fskip) = [];
% 
% % save
% fprintf('Saving to file ...\n');
% save(savestr,'filename','NFRAMES','NCELLS','nv','L','phi','S','calA','meanCalA','stdCalA','calA0','zv','zc','voroAreas','voroCalA');
% 
% 
% fprintf('Thats all! Saved data to %s, ending MATLAB portion.\n',savestr);

end