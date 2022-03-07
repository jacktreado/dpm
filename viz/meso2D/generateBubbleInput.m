function generateBubbleInput(inputfstr, outputfstr)
%% FUNCTION to generate an initial condition with void bubble particles

% read in data
mesoData = readMesoNetwork2D(inputfstr);

% parse for void polygons
NFRAMES = mesoData.NFRAMES;
if NFRAMES ~= 1
    error('generateBubbleInput:tooManyFrames','In generateBubbleInput, inputstr = %s, NFRAMES = %d, which is too many (should be 1). Ending.\n',inputfstr,NFRAMES);
end

% get void polygons of initial condition
polys = voidPolys(mesoData);
NPOLYS = length(polys);

% effective pinned particle sizes and centers
NCELLS = mesoData.NCELLS;
nv = mesoData.nv;
x = mesoData.x;
y = mesoData.y;
r = mesoData.r;
L = mesoData.L;
cx = cellfun(@mean,x)';
cy = cellfun(@mean,y)';
cr = zeros(NCELLS,1);
for nn = 1:NCELLS
    rx = x{nn} - cx(nn);
    ry = y{nn} - cy(nn);
    rads = sqrt(rx.^2 + ry.^2);
    xedge = x{nn} + 3*r{nn}.*(rx./rads);
    yedge = y{nn} + 3*r{nn}.*(ry./rads);
    ex = xedge - cx(nn);
    ey = yedge - cy(nn);
    r2 = sum(ex.^2 + ey.^2)/nv(nn);
    cr(nn) = sqrt(r2);
end
xall = cell2mat(x');
yall = cell2mat(y');
rall = cell2mat(r');
% cx = [cx; xall];
% cy = [cy; yall];
% cr = [cr; rall];

% effective bubble info based on bubble centers
bcx = zeros(NPOLYS,1);
bcy = zeros(NPOLYS,1);
bcr = 0.5*min(rall)*ones(NPOLYS,1);

% loop over polygons, place deformable particles in the center
for pp = 1:NPOLYS
    % void polygon
    polytmp = polys{pp};
    
    % place particle in center of polygon
    pcx = mean(polytmp(:,1));
    pcy = mean(polytmp(:,2));
    
    % use FIRE to place given fixed pinned particles
    fprintf('* Placing bubble %d in polygon center, using FIRE to relax...\n',pp);
    [bcx(pp), bcy(pp)] = tracerFIRE(pcx,pcy,bcr(pp),cx,cy,cr,L);
%     bcx(pp) = pcx;
%     bcy(pp) = pcy;
    fprintf('*...bubble placed! x=%0.4g, y=%0.4f, dx=%0.4g, dy=%0.4g\n',bcx(pp),bcy(pp),pcx-bcx(pp),pcy-bcy(pp));
end

% deformable bubble information
nvb = round(mean(nv))*ones(NPOLYS,1);
xb = cell(NPOLYS,1);
yb = cell(NPOLYS,1);
rb = zeros(NPOLYS,1);
a0b = zeros(NPOLYS,1);
l0b = zeros(NPOLYS,1);
pb = zeros(NPOLYS,1);

% get bubble geometric info
for pp = 1:NPOLYS
    bxtmp = zeros(nvb(pp),1);
    bytmp = zeros(nvb(pp),1);
    brtmp = bcr(pp);
    for vv = 1:nvb(pp)
        bxtmp(vv) = brtmp*cos((2.0*pi*(vv-1))/nvb(pp)) + bcx(pp);
        bytmp(vv) = brtmp*sin((2.0*pi*(vv-1))/nvb(pp)) + bcy(pp);
    end
    xb{pp} = bxtmp;
    yb{pp} = bytmp;
    a0b(pp) = polyarea(bxtmp,bytmp);
    blxtmp = bxtmp([2:nvb(pp) 1]) - bxtmp;
    blytmp = bytmp([2:nvb(pp) 1]) - bytmp;
    bltmp = sqrt(blxtmp.^2 + blytmp.^2);
    pb(pp) = sum(bltmp);
    l0b(pp) = mean(bltmp);
    rb(pp) = 0.5*l0b(pp);
end

% reopen input file for copying
infid = fopen(inputfstr);

% print to file
fid = fopen(outputfstr,'w');
fprintf(fid,'NEWFR\n');
fprintf(fid,'NUMCL %d\n',NCELLS);
fprintf(fid,'NUMBL %d\n',NPOLYS);

% copy input file to output file, except for first two lines
ll = 1;
while ~feof(infid)
    linetmp = fgetl(infid);
    if ll < 3 || contains(linetmp,'ENDFR')
        ll = ll + 1;
        continue;
    end
    fprintf(fid,[linetmp '\n']);
    ll = ll + 1;
end

% now add lines about bubbles
for bb = 1:NPOLYS
    fprintf(fid,'CINFO %6d %8d %20.12g %20.12g %20.12g\n',nvb(bb),0,a0b(bb),a0b(bb),pb(bb));
    for vv = 1:nvb(bb)
        fprintf(fid,'VINFO %6d %8d %20.12g %20.12g %20.12g %20.12g 0 0 0\n',bb-1+NCELLS,vv-1,xb{bb}(vv),yb{bb}(vv),rb(bb),l0b(bb));
    end
end
fprintf(fid,'ENDFR\n');


fclose(infid);
fclose(fid);

end




%% FUNCTION to get void polygons based on data from simulation

function polys = voidPolys(mesoData)

% parse sim info
N = mesoData.NCELLS;
NV = mesoData.nv;
L = mesoData.L;
xtmp = mesoData.x;
ytmp = mesoData.y;
rtmp = mesoData.r;

% get cell-cell contacts
cij = getCellCtcs(N,NV,xtmp,ytmp,rtmp,L);

% call meso void function to get void polygons
cxtmp = cellfun(@mean,xtmp)';
cytmp = cellfun(@mean,ytmp)';
[polys, ~] = getMesoVoidPolygons(cxtmp,cytmp,cij,L(1));

end



%% FUNCTION to determine cell-cell contacts

function cij = getCellCtcs(NCELLS,nv,x,y,r,L)
% initialize output
cij = zeros(NCELLS);

% loop over all vertex-vertex pairs
for ci = 1:NCELLS
    xi = x{ci};
    yi = y{ci};
    ri = r{ci};
    ni = nv(ci);
    for cj = (ci+1):NCELLS
        xj = x{cj};
        yj = y{cj};
        rj = r{cj};
        nj = nv(cj);
        ctcfound = 0;
        for vi = 1:ni
            for vj = 1:nj
                dx = xj(vj) - xi(vi);
                dx = dx - L(1)*round(dx/L(1));
                
                dy = yj(vj) - yi(vi);
                dy = dy - L(2)*round(dy/L(2));
                
                rij = sqrt(dx*dx + dy*dy);
                sij = ri(vi) + rj(vj);
                if (rij < sij)
                    ctcfound = 1;
                    break;
                end
            end
            if ctcfound == 1
                break;
            end
        end
        if ctcfound == 1
            cij(ci,cj) = 1;
            cij(cj,ci) = 1;
        end
    end
end

end





%% FUNCTION to place pin using FIRE

function [x, y] = tracerFIRE(x,y,r,px,py,pr,L)
% pin information
NPINS = length(px);

% force tolerance
Ftol = 1e-12;

% FIRE VARIABLES (hard code in)
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dt          = 0.005;
dtmax       = 10*dt;
dtmin       = 0.02*dt;
NNEGMAX     = 500;
NDELAY      = 20;
npPos       = 0;
npNeg       = 0;
npPMin      = 0;
alpha       = alpha0;

% tracer-pin spring constant
ktp = 0.01;

% force check
fcheck = 10*Ftol;

% initialize forces
fx = 0.0;
fy = 0.0;

% initialize velocities
vx = 0.0;
vy = 0.0;

% USE FIRE to relax forces
it = 0;
itmax = 1e5;

% loop over FIRE protocol while fcheck and kcheck are above minimum
while(fcheck > Ftol && it < itmax)
    % update iterate
    it = it + 1;
    
    % Step 1. calculate P, fnorm, vnorm
    P = sum(fx.*vx) + sum(fy.*vy);

    % plot FIRE information 
    if mod(it,5000) == 0
        fprintf('\nOn FIRE step %d\n',it);
        fprintf('\t ** F = %0.5g\n',fcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** P = %0.5g\n',P);
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** npPMin = %d\n',npPMin);
    end
    
    % Step 2. adjust simulation based on net motion of system
    if P > 0
        % increase positive counter
        npPos = npPos + 1;
        
        % reset negative counter
        npNeg = 0;
        
        % alter simulation if enough positive steps have been taken
        if (npPos > NDELAY)
            % change time step
            if (dt*finc < dtmax)
                dt = dt*finc;
            end
        end
        
        % decrease alpha
        alpha = alpha*falpha;
    else
        % reset positive counter
        npPos = 0;
        
        % increase negative counter
        npNeg = npNeg + 1;
        
        % check for stuck simulation
        if (npNeg > NNEGMAX)
            fprintf('Simulation negative for too long, ending program here.\n')
            error('FIRE did not converge');
        end
        
        % decrease time step if past initial delay
        if (it > NDELAY)
            % decrease time step
            if (dt*fdec > dtmin)
                dt = dt*fdec;
            end
            
            % reset alpha
            alpha = alpha0;
        end
        
        % take a half step backwards
        x = x - 0.5*dt*vx;
        y = y - 0.5*dt*vy;
        
        % reset velocities to 0
        vx = 0.0;
        vy = 0.0;
    end
    
    % Step 1 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;
    
    % update vnorm and fnorm for FIRE
    vnorm = sqrt(sum(vx.*vx) + sum(vy.*vy));
    fnorm = sqrt(sum(fx.*fx) + sum(fy.*fy));

    
    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*(fx./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*(fy./fnorm)*vnorm;
    end
    
    % do first verlet update for vertices (assume unit mass)
    x = x + dt*vx;
    y = y + dt*vy;
    
    % compute forces on pin
    fx = 0.0;
    fy = 0.0;
    for pp = 1:NPINS
        dx = px(pp) - x;
        dx = dx - L(1)*round(dx/L(1));
        dy = py(pp) - y;
        dy = dy - L(2)*round(dy/L(2));
        rij = sqrt(dx*dx + dy*dy);
        sij = pr(pp) + r;
        if rij < sij
            ftmp = (1 - (rij/sij))/sij;
            ux = dx/rij;
            uy = dy/rij;
            fxtmp = ftmp*ux;
            fytmp = ftmp*uy;
            fx = fx - fxtmp;
            fy = fy - fytmp;
        end
    end
    
    % step 3 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % update force check
    fcheck = sqrt(sum(fx.^2 + fy.^2));
end

end