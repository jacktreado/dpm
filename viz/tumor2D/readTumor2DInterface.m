
function tumorConfigData = readTumor2DInterface(fstr)
%% FUNCTION to read in DPM config data given file string
% Specifically for tumor interface invasion data
% NOTE: ASSUMING FIXED NUMBER OF CELLS


% print info to console
finfo = dir(fstr);
fprintf('-- Reading in %s\n',finfo.name);
fprintf('-- File size = %f MB\n',finfo.bytes/1e6);

% open file stream
fid = fopen(fstr);

% read in sim details from first frame
ndata       = textscan(fid,'NUMCL %f %f',1,'HeaderLines',1);
NCELLS      = ndata{1};
tN          = ndata{2};
fline       = fgetl(fid);


dttmp       = textscan(fid,'TSTEP %f',1);
fline       = fgetl(fid);

ttmp        = textscan(fid,'TCURR %f',1);
fline       = fgetl(fid);

phitmp      = textscan(fid,'PACKF %f',1);
fline       = fgetl(fid);

Ltmp        = textscan(fid,'BOXSZ %f %f',1);
fline       = fgetl(fid);

Stmp        = textscan(fid,'STRSS %f %f %f',1);
fline       = fgetl(fid);

WPtmp       = textscan(fid,'WPRSS %f %f %f',1);

% cells to save 
NFRAMES = 1e6;

phi     = zeros(NFRAMES,1);
dt      = zeros(NFRAMES,1);
t       = zeros(NFRAMES,1);
L       = zeros(NFRAMES,2);
S       = zeros(NFRAMES,3);
WP      = zeros(NFRAMES,2);

nv      = zeros(NFRAMES,NCELLS);
zc      = zeros(NFRAMES,NCELLS);
zv      = zeros(NFRAMES,NCELLS);
a0      = zeros(NFRAMES,NCELLS);
a       = zeros(NFRAMES,NCELLS);
p       = zeros(NFRAMES,NCELLS);

x       = cell(NFRAMES,NCELLS);
y       = cell(NFRAMES,NCELLS);
r       = cell(NFRAMES,NCELLS);
l0      = cell(NFRAMES,NCELLS);
t0      = cell(NFRAMES,NCELLS);


% tumor-specific
psi     = zeros(NFRAMES,tN);
Dr      = zeros(NFRAMES,tN);

% adipocyte-specific
px      = zeros(NFRAMES,NCELLS-tN);
py      = zeros(NFRAMES,NCELLS-tN);
attch   = false(NFRAMES,NCELLS-tN);

% number of frames found
nf = 1;

% loop over frames, read in data
while ~feof(fid)
    % get time info
    dt(nf) = dttmp{1};
    t(nf) = ttmp{1};
    
    % get packing fraction
    phi(nf) = phitmp{1};
    
    % get box length
    L(nf,1) = Ltmp{1};
    L(nf,2) = Ltmp{2};
    
    % get stress info
    S(nf,1) = Stmp{1};
    S(nf,2) = Stmp{2};
    S(nf,3) = Stmp{3}; 
    
    % get wall pressure
    WP(nf,1) = WPtmp{1};
    WP(nf,2) = WPtmp{2};
    
    % get info about deformable particle
    for nn = 1:NCELLS
        % get cell pos and asphericity
        if nn <= tN
            cInfoTmp = textscan(fid,'CINFO %f %f %f %f %f %f %f %f',1);   
            fline = fgetl(fid);     % goes to next line in file
            
            psi(nf,nn) = cInfoTmp{7};
            Dr(nf,nn) = cInfoTmp{8};
        else
            cInfoTmp = textscan(fid,'CINFO %f %f %f %f %f %f %f %f %f',1);   
            fline = fgetl(fid);     % goes to next line in file
            
            px(nf,nn-tN) = cInfoTmp{7};
            py(nf,nn-tN) = cInfoTmp{8};
            attch(nf,nn-tN) = cInfoTmp{9};
        end
        NVTMP = cInfoTmp{1};
        nv(nf,nn) = NVTMP;
        zc(nf,nn) = cInfoTmp{2};
        zv(nf,nn) = cInfoTmp{3};
        a0(nf,nn) = cInfoTmp{4};
        a(nf,nn) = cInfoTmp{5};
        p(nf,nn) = cInfoTmp{6};
        
        % get vertex positions
        vInfoTmp = textscan(fid,'VINFO %*f %*f %f %f %f %f %f',NVTMP); 
        fline = fgetl(fid);     % goes to next line in file

        % parse data
        x{nf,nn} = vInfoTmp{1};
        y{nf,nn} = vInfoTmp{2};
        r{nf,nn} = vInfoTmp{3};
        l0{nf,nn} = vInfoTmp{4};
        t0{nf,nn} = vInfoTmp{5};
    end
    
    % increment frame count
    nf = nf + 1;
    
    % get ENDFR string, check that read is correct
    fline = fgetl(fid);
    endFrameStr = sscanf(fline,'%s');
    
    if ~strcmp(endFrameStr,'ENDFR')
        error('Miscounted number of particles in .pos file, ending data read in');
    end
    
    % start new frame information
    if ~feof(fid)
        % get NEWFR line
        fline       = fgetl(fid);
        newFrameStr = sscanf(fline,'%s %*s');
        if ~strcmp(newFrameStr,'NEWFR')
            error('NEWFR not encountered when expected when heading to new frame...check line counting. ending.');
        end
        
        % read in sim details from first frame
        NCELLStmp       = textscan(fid,'NUMCL %f %f');
        
        % read in next frame header data
        dttmp       = textscan(fid,'TSTEP %f',1);
        fline       = fgetl(fid);

        ttmp        = textscan(fid,'TCURR %f',1);
        fline       = fgetl(fid);

        phitmp      = textscan(fid,'PACKF %f',1);
        fline       = fgetl(fid);

        Ltmp        = textscan(fid,'BOXSZ %f %f',1);
        fline       = fgetl(fid);

        Stmp        = textscan(fid,'STRSS %f %f %f',1);
        fline       = fgetl(fid);

        WPtmp       = textscan(fid,'WPRSS %f %f %f',1);
        fline       = fgetl(fid);
    end
end

% delete extra information
if (nf < NFRAMES)
    NFRAMES = nf-1;
    phi(nf:end) = [];
    dt(nf:end) = [];
    t(nf:end) = [];
    L(nf:end,:) = [];
    S(nf:end,:) = [];
    WP(nf:end,:) = [];
    
    nv(nf:end,:) = [];
    zc(nf:end,:) = [];
    zv(nf:end,:) = [];
    a0(nf:end,:) = [];
    a(nf:end,:) = [];
    p(nf:end,:) = [];
    psi(nf:end,:) = [];
    Dr(nf:end,:) = [];
    px(nf:end,:) = [];
    py(nf:end,:) = [];
    attch(nf:end,:) = [];
    
    x(nf:end,:) = [];
    y(nf:end,:) = [];
    r(nf:end,:) = [];
    l0(nf:end,:) = [];
    t0(nf:end,:) = [];
end

% close position file
fclose(fid);

% store cell pos data into struct
tumorConfigData                 = struct('NFRAMES',NFRAMES,'NCELLS',NCELLS,'tN',tN);
tumorConfigData.phi             = phi;
tumorConfigData.dt              = dt;
tumorConfigData.t               = t;
tumorConfigData.L               = L;
tumorConfigData.S               = S;
tumorConfigData.WP              = WP;

tumorConfigData.nv              = nv;
tumorConfigData.zc              = zc;
tumorConfigData.zv              = zv;
tumorConfigData.a0              = a0;
tumorConfigData.a               = a;
tumorConfigData.p               = p;
tumorConfigData.psi             = psi;
tumorConfigData.Dr              = Dr;
tumorConfigData.px              = px;
tumorConfigData.py              = py;
tumorConfigData.attch           = attch;

tumorConfigData.x               = x;
tumorConfigData.y               = y;
tumorConfigData.r               = r;
tumorConfigData.l0              = l0;
tumorConfigData.t0              = t0;

end