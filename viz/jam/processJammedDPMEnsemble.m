function processJammedDPMEnsemble(ensembleStr, savestr)
%% Process ensemble of jammed dpm configurations
% ** ASSUMING ONLY 1 FRAME, will skip configs that are empty or have
% multiple frames

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
filename    = cell(NEN,1);
phi         = zeros(NEN,1);
S           = zeros(NEN,3);
calA        = cell(NEN,1);
meanCalA    = zeros(NEN,1);
stdCalA     = zeros(NEN,1);
calA0       = cell(NEN,1);
zc          = cell(NEN,1);
zv          = cell(NEN,1);

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
    filename{ee} = filename;
    
    % load in data
    dpmConfigData = readDPMConfig(fstr);
    
    % check # of frames
    NFRAMES = dpmConfigData.NFRAMES;
    if NFRAMES ~= 1
        fprintf('** File %s has %d frames, skipping...\n',fname,NFRAMES);
        fskip(ee) = true;
        continue;
    end
    NCELLS = dpmConfigData.NCELLS;
    
    % save data
    phi(ee) = dpmConfigData.phi;
    
    Stmp    = dpmConfigData.S;
    S(ee,1) = Stmp(1);
    S(ee,2) = Stmp(2);
    S(ee,3) = Stmp(3);
    
    ptmp = dpmConfigData.p;
    atmp = dpmConfigData.a;
    calA{ee} = ptmp.^2./(4.0*pi*atmp);
    meanCalA(ee) = mean(calA{ee});
    stdCalA(ee) = std(calA{ee});
    
    l0tmp = dpmConfigData.l0;
    a0tmp = dpmConfigData.a0;
    calA0tmp = zeros(NCELLS,1);
    for nn = 1:NCELLS
        p0tmp = sum(l0tmp{nn});
        calA0tmp(nn) = p0tmp^2/(4.0*pi*a0tmp(nn));
    end
    calA0{ee} = calA0tmp;
    
    zc{ee} = dpmConfigData.zc;
    zv{ee} = dpmConfigData.zv;
end

% delete extra entries
filename(fskip) = [];
phi(fskip) = [];
S(fskip,:) = [];
calA(fskip) = [];
meanCalA(fskip) = [];
stdCalA(fskip) = [];
calA0(fskip) = [];
zv(fskip) = [];
zc(fskip) = [];


% save
save(savestr,'filename','phi','S','calA','meanCalA','stdCalA','calA0','zv','zc');

end