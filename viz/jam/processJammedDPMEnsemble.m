function processJammedDPMEnsemble(ensembleLoc, savestr)
%% Process ensemble of jammed dpm configurations
% ** ASSUMING ONLY 1 FRAME, will skip configs that are empty or have
% multiple frames

% get files that make up ensemble
if strcmp(ensembleLoc(end),'/')
    ensembleLoc(end) = [];
end
ensList = flist([ensembleLoc '/*.pos']);
NEN = length(ensList);
if NEN == 0
    error('processJammedDPMEnsemble:noFilesFound','\nNo files found in ensemble location %s. Ending here.\n',ensembleLoc);
else
    fprintf('** Processing %d files stored in %s\n',NEN,ensembleLoc);
end

%% Loop over files in ensemble, extract data

% files to skip
fskip       = false(NEN,1);

% data to save
phi         = zeros(NEN,1);
S           = zeros(NEN,3);
calA        = cell(NEN,1);
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
    
    % load in data
    dpmConfigData = readDPMConfig(fstr);
    
    % check # of frames
    NFRAMES = dpmConfigData.NFRAMES;
    if NFRAMES ~= 1
        fprintf('** File %s has %d frames, skipping...\n',fname,NFRAMES);
        fskip(ee) = true;
        continue;
    end
    
    % save data
    phi(ee) = dpmConfigData.phi;
    
    Stmp    = dpmConfigData.S;
    S(ee,1) = Stmp(1);
    S(ee,2) = Stmp(2);
    S(ee,3) = Stmp(3);
    
    ptmp = dpmConfigData.p;
    atmp = dpmConfigData.a;
    calA{ee} = ptmp.^2./(4.0*pi*atmp);
    
    l0tmp = dpmConfigData.l0;
    
end


end