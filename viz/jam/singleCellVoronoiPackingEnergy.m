function U = singleCellVoronoiPackingEnergy(nv,x,y,r,a0,l0,V)
%% FUNCTION to compute energy from DPM particle packed into svoro cell

% r-squared (for quick overlap check)
r2 = r.^2;

% voronoi coordinates and geometric info
nsvv = size(V,1);
svx = V(:,1);
svy = V(:,2);

% compute shape parameter of voronoi
svlx = svx([2:nsvv 1]) - svx;
svly = svy([2:nsvv 1]) - svy;
svl = sqrt(svlx.^2 + svly.^2);
svp = sum(svl);

% get cell edge unit vectors
usvlx = svlx./svl;
usvly = svly./svl;

% get cell normals
svnx = -usvly;
svny = usvlx;

%% FIRE relaxation

% fixed time parameters
dt0         = 0.005;
itmax       = 1e7;
plotskip    = 5000;
Ftol        = 1e-12;

% indexing
im1 = [nv 1:nv-1];
ip1 = [2:nv 1];

% force parameters
rho0            = sqrt(a0);                 % units: length
fa              = 1/rho0;                   % units: inv length, because of grad a
fl              = rho0/l0;                  % units: dim. less

% FIRE VARIABLES (hard code in)
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dtmax       = 10*dt0;
dtmin       = 0.02*dt0;
dt          = dt0;
NNEGMAX     = 500;
NDELAY      = 20;
npPos       = 0;
npNeg       = 0;
npPMin      = 0;
alpha       = alpha0;

% force check
fcheck = 10*Ftol;

% initialize forces
fx = zeros(nv,1);
fy = zeros(nv,1);

% initialize velocities
vx = zeros(nv,1);
vy = zeros(nv,1);

% energy
U = 0;
Ua = 0;
Ul = 0;
Usv = 0;

% USE FIRE to relax forces
it = 0;

% loop over FIRE protocol while fcheck and kcheck are above minimum
while(fcheck > Ftol && it < itmax)
    % update iterate
    it = it + 1;
    
    % Step 1. calculate P, fnorm, vnorm
    P = sum(fx.*vx) + sum(fy.*vy);

    % plot FIRE information 
    if mod(it,plotskip) == 0
        fprintf('& ');
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
        vx = zeros(nv,1);
        vy = zeros(nv,1);
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
    
    % * * * * * * * * * * * * * * * * * *
    % calculate forces based on positions
    % * * * * * * * * * * * * * * * * * *
    
    % reset forces
    fx = zeros(nv,1);
    fy = zeros(nv,1);
    
    % -- perimeter force
    
    % segment vectors
    lvx = x(ip1) - x;
    lvy = y(ip1) - y;
    
    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    
    % segment unit vectos
    ulvx = lvx./l;
    ulvy = lvy./l;
    
    % segment strain
    dli = (l./l0) - 1.0;
    dlim1 = (l(im1)./l0) - 1.0;
    
    % perimeter force at this iteration
    flx = fl*(dli.*ulvx - dlim1.*ulvx(im1));
    fly = fl*(dli.*ulvy - dlim1.*ulvy(im1));
    
    % add to total force
    fx = fx + flx;
    fy = fy + fly;
    
    % add to potential energy
    Ul = 0.5*sum(dli.^2);
    
    % -- area force
    a = polyarea(x, y);
    areaStrain = (a/a0) - 1.0;
    
    % area force at this iteration
    fax = fa*0.5*areaStrain.*(y(im1) - y(ip1));
    fay = fa*0.5*areaStrain.*(x(ip1) - x(im1));
    
    % add to total force
    fx = fx + fax;
    fy = fy + fay;
    
    % add to potential energy
    Ua = 0.5*areaStrain^2;
    
    
    % -- force due to walls on each vertex
    
    % distances to vertices
    dvx = x' - svx;
    dvy = y' - svy;
    dv2 = dvx.^2 + dvy.^2;
    dv = sqrt(dv2);
    
    % unit vectors to vertices
    udvx = dvx./dv;
    udvy = dvy./dv;
    
    % get sstar (closest point to entire line)
    sstar = (dv./svl).*(udvx.*usvlx + udvy.*usvly);
    sedge = sstar;
    
    % for sstar outside 0 and 1, check closes point to edge
    offBackEdge = sstar < 0;
    offFrontEdge = sstar > 1;
    sedge(offBackEdge(:)) = zeros(sum(offBackEdge(:)),1);
    sedge(offFrontEdge(:)) = ones(sum(offFrontEdge(:)),1);
    
    % get distances to closest point on edge
    h2 = dv2 + (svl.^2).*sedge.*(sedge - 2.*sstar);

    % overlapping indices
    ovinds = (h2 < r2');

    % add to forces
    Usv = 0;
    for nn = 1:nv
        ovtmp = ovinds(:,nn);
        h = sqrt(h2(ovtmp,nn));
        rtmp = r(nn);
        nov = sum(ovtmp);
        if (nov > 0)
            % add to forces
            ftmp = (1 - (h./rtmp))./(rtmp*svp*nov);
            fedgex = ftmp.*svnx(ovtmp);
            fedgey = ftmp.*svny(ovtmp);
            fx(nn) = fx(nn) + sum(fedgex);
            fy(nn) = fy(nn) + sum(fedgey);
            
            % add to energy
            Usv = Usv + (0.5/(svp*nov))*sum((1 - (h./rtmp)).^2);
        end
    end
    
    % step 3 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % update force check
    fcheck = sqrt(sum(fx.^2 + fy.^2)/nv);
    U = Ua + Usv + Ul;
end
if (it == itmax)
    U = -1;
    fprintf('DID NOT FIND MINIMUM, setting final U = %0.4g\n',U);
else
    fprintf('U=%0.4g  ',U);
end

end