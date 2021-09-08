function processMesoPenta2D(floc,fpattern,savestr)
%% FUNCTION to analyze mesophyll network simulations

% get list of files that fit pattern
flist = dir([floc '/' fpattern '*.pos']);
NSIMS = length(flist);
if NSIMS == 0
    error('processMesoNetwork2D:noFilesFound','No files found with pattern %s in location %s\n',fpattern,floc);
else
    fprintf('Found %d file using pattern %s, processing...\n',NSIMS,fpattern);
end

% number of cells (constant)
NCELLS = 6;


%% Loop over files, save information

% files to skip
fskip = false(NSIMS,1);

% data to save
params = zeros(NSIMS,10);
netx = cell(NSIMS,1);
nety = cell(NSIMS,1);
zcList = cell(NSIMS,1);
calAList = cell(NSIMS,1);
calA0List = cell(NSIMS,1);

% Loop over sims, extract parameters + configurations
for ss = 1:NSIMS
	% get file info
	fname = flist(ss).name;
	floc = flist(ss).folder;
	fsize = flist(ss).bytes;
	if fsize == 0
		fprintf('** File %s has no data, skipping...\n',fname);
		fskip(ss) = true;
	else
		fprintf('** Reading in data from file %s ... \n',fname);

	end
	fstr = [floc '/' fname];
	pentaData = readMesoPin2D(fstr);

	% extract parameters from file name
    paramtmp = sscanf(fname,'penta_n%f_ca%f_kb0%f_be%f_cd%f_ch%f_cL%f_aL%f_cB%f_cKb%f.pos');
    params(ss,:) = paramtmp';
    
	% parse data from file
    NFRAMES = pentaData.NFRAMES;
	x = pentaData.x;
    y = pentaData.y;
    zc = pentaData.zc;
    a0 = pentaData.a0;
    l0 = pentaData.l0;
    p = pentaData.p;
    a = pentaData.a;
    
    % get connection
    zcList{ss} = zc;

    % get shapes
    calA0 = zeros(NFRAMES,NCELLS);
    calA = zeros(NFRAMES,NCELLS);
    for ff = 1:NFRAMES
        a0tmp = a0(ff,:);
        l0tmp = l0(ff,:);
        for cc = 1:NCELLS
            p0tmp = sum(l0tmp{cc});
            calA0(ff,cc) = p0tmp^2/(4.0 * pi * a0tmp(cc));
            calA(ff,cc) = p(ff,cc)^2/(4.0 * pi * a(ff,cc));
        end
    end
    calA0List{ss} = calA0;
    calAList{ss} = calA;
    
    
    % get last network config
    idx = zc == 3;
    cnxidx = find(sum(idx(:,2:end),2) == 5);
    lastcnxidx = cnxidx(end);
    netx{ss} = x(lastcnxidx,:);
    nety{ss} = y(lastcnxidx,:);
end

% delete extra entries, save
flist(fskip) = [];
params(fskip,:) = [];
zcList(fskip) = [];
calA0(fskip) = [];
calA(fskip) = [];
netx(fskip) = [];
nety(fskip) = [];

% save
save(savestr,'flist','params','zcList','calA0','calA','netx','nety');
end