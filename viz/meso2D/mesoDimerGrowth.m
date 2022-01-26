function mesoDimerGrowth(savestr,NGROWTH,kl,kb,kc,del_l0,P0)
%% FUNCTION to grow meso cell dimer in box with various parameters

% simulation details
n           = 32;
a0base      = 1.0;
calA0       = n*tan(pi/n)/pi;
l0base      = sqrt(4.0*pi*calA0*a0base)/n;
th0base     = 0.0;

% vertex-based shape parameters
r           = (0.5*l0base)*ones(n,2);
a0          = a0base*ones(1,2);
l0          = l0base*ones(n,2);
th0         = th0base*ones(n,2);

% arrays for degrees of freedom
x           = zeros(n,2);
y           = zeros(n,2);

% compression parameter (0 = no compression, 1 = total compression)
cmp         = 0.2;
dL          = 1.0 - cmp;


%% Simulation setup

% position of first particle (top-half of box)
dth = (2.0*pi)/n;
th  = 0:dth:(2.0*pi - dth);
x0  = cos(th)';
y0	= sin(th)';
ainit = polyarea(x0,y0);
x0 = x0./sqrt(ainit);
y0 = y0./sqrt(ainit);
r0 = sqrt(x0(1)^2 + y0(1)^2);

% maximum box lengths
Lxmax = 2.0*r0 + l0base;
Lymax = 2.0*Lxmax;

% set positions of vertices
x(:,1) = x0 + 0.5*Lxmax;
y(:,1) = y0 + 3*r0 + 1.5*l0base;

x(:,2) = x0 + 0.5*Lxmax;
y(:,2) = y0 + r0 + 0.5*l0base;

% set box length based on compression
Lx = Lxmax;
Ly = dL*Lymax;
L = zeros(1,2);
L(1) = Lx;
L(2) = Ly;

% strain positions
x = x.*(Lx./Lxmax);
y = y.*(Ly./Lymax);

% initialize bonded matrix
bij = zeros(n);
bw = zeros(n,3,2);


%% Relax at constant volume for initial condition

% find relaxed configuration 
[x,y,fx,fy] = mesoDimerEMin(x,y,r,a0,l0,th0,L,kl,kb,kc,kc,bij,bw);

% measure pressure
Pext = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kc,bij,bw);
Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
P = Pext + Pint;
fprintf('Initial pressure without bonds at compression cmp = %0.5g;  P = %0.5g\n',cmp,P);

%% Setup bond network

% initialize pairwise bonds (+1)
for ii = 1:n
    xi = x(ii,1);
    yi = y(ii,1);
    ri = r(ii,1);
    for jj = 1:n
        xj = x(jj,2);
        yj = y(jj,2);
        rj = r(jj,2);
        sij = ri + rj;
        dx = xj - xi;
        dy = yj - yi;
        rij = sqrt(dx*dx + dy*dy);
        if rij < sij
            bij(ii,jj) = 1;
        end
    end
end

% initialize wall bonds (-1)
for nn = 1:2
    for ii = 1:n
        % left bond = 1
        if x(ii,nn) < r(ii,nn)
            bw(ii,1,nn) = 1;
            bw(ii,2,nn) = 0;
            bw(ii,3,nn) = y(ii,nn);
        end
        % bottom bond = 2
        if y(ii,nn) < r(ii,nn)
            bw(ii,1,nn) = 2;
            bw(ii,2,nn) = x(ii,nn);
            bw(ii,3,nn) = 0;
        end
        % right bond = 3
        if x(ii,nn) > L(1) - r(ii,nn)
            bw(ii,1,nn) = 3;
            bw(ii,2,nn) = L(1);
            bw(ii,3,nn) = y(ii,nn);
        end
        % top bond = 4
        if y(ii,nn) > L(2) - r(ii,nn)
            bw(ii,1,nn) = 4;
            bw(ii,2,nn) = x(ii,nn);
            bw(ii,3,nn) = L(2);
        end
    end
end

%% Grow + minimize enthalpy

% initial enthalpy min
fprintf('** Initial enthalpy minimization...\n');
[x,y,fx,fy,L,bw] = mesoDimerHMin(x,y,r,a0,l0,th0,P0,L,kl,kb,kc,kc,bij,bw);
Pext = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kc,bij,bw);
Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
P = Pext + Pint;
fprintf('...initial pressure P = %0.5g, P0 = %0.5g\n',P,P0);

% save during growth steps
xList = cell(NGROWTH,2);
yList = cell(NGROWTH,2);
LList = zeros(NGROWTH,2);
a0List = zeros(NGROWTH,2);
l0List = cell(NGROWTH,2);
calAList = zeros(NGROWTH,2);
calA0List = zeros(NGROWTH,2);

% Loop over growth and relaxation steps
for gg = 1:NGROWTH
    % grow void perimeter
    for nn = 1:2
        for ii = 1:n
            if bw(ii,1,nn) == 0 && ((nn == 1 && sum(bij(ii,:)) == 0) || (nn == 2 && sum(bij(:,ii)) == 0))
                l0(ii,nn) = l0(ii,nn)*(1 + del_l0);
            end
        end
    end
                
    % enthalpy minimize
    [x,y,fx,fy,L,bw] = mesoDimerHMin(x,y,r,a0,l0,th0,P0,L,kl,kb,kc,kc,bij,bw);
    
    % compute pressure, print
    Pext = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kc,bij,bw);
    Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
    P = Pext + Pint;
    lx = x([2:n 1],:) - x;
    ly = y([2:n 1],:) - y;
    l = sqrt(lx.^2 + ly.^2);
    a = polyarea(x,y);
    da = (a./a0) - 1;
    p = sum(l,1);
    calAtmp = p.^2./(4.0*pi*a);
    calA0tmp = sum(l0).^2./(4.0*pi*a0); 
    
    fprintf('\n** Enthalpy minimized + printing\n');
    fprintf('** Pext = %0.5g\n',Pext);
    fprintf('** Pint = %0.5g\n',Pint);
    fprintf('** P = %0.5g\n',P);
    fprintf('** L = %0.5g\n',L(1));
    fprintf('** calA1 = %0.5g, calA2 = %0.5g\n',calAtmp(1),calAtmp(2));
    fprintf('** calA01 = %0.5g, calA02 = %0.5g\n',calA0tmp(1),calA0tmp(2));
    fprintf('** phi = %0.5g\n',sum(a)/(L(1)*L(2)));
    fprintf('** area strain = %0.5g, %0.5g\n\n\n',da(1),da(2));
    
    % save relaxed data
    xList{gg,1} = x(:,1);
    xList{gg,2} = x(:,2);
    yList{gg,1} = y(:,1);
    yList{gg,2} = y(:,2);
    LList(gg,:) = L;
    a0List(gg,:) = a0;
    l0List{gg,1} = l0(:,1);
    l0List{gg,2} = l0(:,2);
    calAList(gg,:) = calAtmp;
    calA0List(gg,:) = calA0tmp;
    
    % set preferred area to be instantaneous
    a = polyarea(x,y);
    a0 = a;
end

% save
save(savestr,...
    'xList','yList','LList','a0List',...
    'l0List','calAList','calA0List');

end




%% -- AUXILLIARY FUNCTIONS

function [x,y,fx,fy] = mesoDimerEMin(x,y,r,a0,l0,th0,L,kl,kb,kc,kw,bij,bw)
%% FUNCTION to relax total potential energy using fire for meso dimer

% number of vertices
n = size(x,1);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];

% forces and velocities
vx = zeros(n,2);
vy = zeros(n,2);
fx = zeros(n,2);
fy = zeros(n,2);

% FIRE parameters
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dt0         = 0.01;
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
Ftol = 1e-12;
fcheck = 10*Ftol;

% USE FIRE to relax forces
it = 0;
itmax = 1e4;

% loop over FIRE protocol while fcheck and kcheck are above minimum
while(fcheck > Ftol && it < itmax)
    % update iterate
    it = it + 1;

    % Step 1. calculate P, fnorm, vnorm
    P = sum(fx.*vx,'all') + sum(fy.*vy,'all');

    % plot FIRE information 
    if mod(it,20000) == 0
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
        vx = zeros(n,2);
        vy = zeros(n,2);
    end

    % Step 1 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % update vnorm and fnorm for FIRE
    vnorm = sqrt(sum(vx.*vx,'all') + sum(vy.*vy,'all'));
    fnorm = sqrt(sum(fx.*fx,'all') + sum(fy.*fy,'all'));


    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*(fx./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*(fy./fnorm)*vnorm;
    end

    % do first verlet update for vertices (assume unit mass)
    x = x + dt*vx;
    y = y + dt*vy;
    
    % figure out forces on all vertices for both cells
    fx = zeros(n,2);
    fy = zeros(n,2);
    for nn = 1:2
        % shape parameters
        xtmp = x(:,nn);
        ytmp = y(:,nn);
        rtmp = r(:,nn);
        a0tmp = a0(nn);
        l0tmp = l0(:,nn);
        th0tmp = th0(:,nn);

        % -- perimeter force

        % segment vectors
        lvx = xtmp(ip1) - xtmp;
        lvy = ytmp(ip1) - ytmp;

        % update perimeter segment lengths
        l = sqrt(lvx.^2 + lvy.^2);

        % segment unit vectos
        ulvx = lvx./l;
        ulvy = lvy./l;

        % segment strain
        dli = ((l./l0tmp) - 1.0)./l0tmp;
        dlim1 = ((l(im1)./l0tmp(im1)) - 1.0)./l0tmp(im1);

        % perimeter force at this iteration
        flx = kl*(dli.*ulvx - dlim1.*ulvx(im1));
        fly = kl*(dli.*ulvy - dlim1.*ulvy(im1));

        % add to total force, potential
        fx(:,nn) = fx(:,nn) + flx;
        fy(:,nn) = fy(:,nn) + fly;

        % -- area force
        % update area
        atmp = polyarea(xtmp, ytmp);
        da = (atmp/a0tmp) - 1.0;

        % area force at this iteration
        fax = 0.5*da.*(ytmp(im1) - ytmp(ip1))./a0tmp;
        fay = 0.5*da.*(xtmp(ip1) - xtmp(im1))./a0tmp;

        % add to total force, potential
        fx(:,nn) = fx(:,nn) + fax;
        fy(:,nn) = fy(:,nn) + fay;
        
        % -- bending force

        % get sine + cosine
        si = lvx.*lvy(im1) - lvy.*lvx(im1);
        ci = lvx.*lvx(im1) + lvy.*lvy(im1);
        thi = atan2(si,ci);
        dth = thi - th0tmp;

        % get normal vector
        nix = lvy;
        niy = -lvx;

        % construct force components
        fbix = (dth - dth(ip1)).*nix./(l.^2);
        fbiy = (dth - dth(ip1)).*niy./(l.^2);

        fbim1x = -fbix(im1);
        fbim1y = -fbiy(im1);

        % add up forces
        fbx = kb*(fbim1x + fbix);
        fby = kb*(fbim1y + fbiy);

        % add to force, potential
        fx(:,nn) = fx(:,nn) + fbx;
        fy(:,nn) = fy(:,nn) + fby;
        
        
        % check wall forces
        fx_wall = zeros(n,1);
        fy_wall = zeros(n,1);
        
        % wall indexes
        left_idx = xtmp < rtmp | bw(:,1,nn) == 1;
        bottom_idx = ytmp < rtmp | bw(:,1,nn) == 2;
        right_idx = xtmp > L(1) - rtmp | bw(:,1,nn) == 3;
        top_idx = ytmp > L(2) - rtmp | bw(:,1,nn) == 4;
                
        % get wall forces
        fx_wall(left_idx) = (kw./rtmp(left_idx)).*(1.0 - (xtmp(left_idx)./rtmp(left_idx)));
        fy_wall(bottom_idx) = (kw./rtmp(bottom_idx)).*(1.0 - (ytmp(bottom_idx)./rtmp(bottom_idx)));
        
        fx_wall(right_idx) = -(kw./rtmp(right_idx)).*(1.0 - ((L(1) - xtmp(right_idx))./rtmp(right_idx)));
        fy_wall(top_idx) = -(kw./rtmp(top_idx)).*(1.0 - ((L(2) - ytmp(top_idx))./rtmp(top_idx)));
        
        fx(:,nn) = fx(:,nn) + fx_wall;
        fy(:,nn) = fy(:,nn) + fy_wall;
        
        
        % check contact forces
        if nn == 1
            for ii = 1:n
                xi = xtmp(ii);
                yi = ytmp(ii);
                for jj = 1:n
                    % contact distance
                    sij = rtmp(ii) + r(jj,2);
                    
                    % distance
                    xj = x(jj,2);
                    yj = y(jj,2);
                    
                    dx = xj - xi;
                    dy = yj - yi;
                    
                    rij = sqrt(dx*dx + dy*dy);
                    
                    % check contacts
                    if (rij < sij || bij(ii,jj) == 1)
                        ftmp = (kc/sij)*(1 - (rij/sij));
                        fxtmp = ftmp*(dx/rij);
                        fytmp = ftmp*(dy/rij);
                        
                        fx(ii,1) = fx(ii,1) - fxtmp;
                        fy(ii,1) = fy(ii,1) - fytmp;
                        
                        fx(jj,2) = fx(jj,2) + fxtmp;
                        fy(jj,2) = fy(jj,2) + fytmp;
                    end
                end
            end
        end
    end
    

    % step 3 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % update force check
    fcheck = sqrt(sum(fx.^2 + fy.^2,'all')/(2.0*n));
end



end


function [x,y,fx,fy,L,bw] = mesoDimerHMin(x,y,r,a0,l0,th0,P0,L,kl,kb,kc,kw,bij,bw)
%% FUNCTION to relax total potential energy using fire for meso dimer

% number of vertices
n = size(x,1);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];

% forces and velocities
vx = zeros(n,2);
vy = zeros(n,2);
fx = zeros(n,2);
fy = zeros(n,2);

% enthalpy variables
V = L(1)*L(2);
Pi = 0.0;
P = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
boxAR = L(2)/L(1);
PM = 1;

% FIRE parameters
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dt0         = 0.001;
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
Ftol = 1e-12;
fcheck = 10*Ftol;

% USE FIRE to relax forces
it = 0;
itmax = 1e8;

% loop over FIRE protocol while fcheck and kcheck are above minimum
while(fcheck > Ftol && it < itmax)
    % update iterate
    it = it + 1;
    
    % Step 1 (VV)
    vx = vx + 0.5*dt*(fx - 0.5*vx*(Pi/V));
    vy = vy + 0.5*dt*(fy - 0.5*vy*(Pi/V));
    Pi = Pi + 0.5*dt*(P - P0)/PM;

    % Step 1. calculate P, fnorm, vnorm
    PFIRE = sum(vx.*(fx - 0.5*vx*(Pi/V)),'all') + sum(vy.*(fy - 0.5*vy*(Pi/V)),'all');
    PFIRE = PFIRE + Pi*(P - P0);

    % plot FIRE information 
    if mod(it,20000) == 0
        fprintf('\n\t && On FIRE HMin step %d\n',it);
        fprintf('\t ** F = %0.5g\n',fcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** PFIRE = %0.5g\n',PFIRE);
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** npPMin = %d\n',npPMin);
        fprintf('\t ** V = %0.5g\n',V);
        fprintf('\t ** P = %0.5g\n',P);
        fprintf('\t ** Pi = %0.5g\n',Pi);
        
        %{  
            &&&&&&&&&&&&&&&&&&&&&&
            Plotting for debugging
            &&&&&&&&&&&&&&&&&&&&&&
        
            figure(1), clf, hold on, box on;
            for nn = 1:2
                for vv = 1:n
                    rectangle('Position',[x(vv,nn)-0.5*l0(vv,nn),y(vv,nn)-0.5*l0(vv,nn),l0(vv,nn),l0(vv,nn)],'Curvature',[1 1],'FaceColor','b');
                end
            end
            axis equal;
            ax = gca;
            ax.XTick = [];
            ax.YTick = [];
            ax.XLim = [-0.25 1.25]*L(1);
            ax.YLim = [-0.25 1.25]*L(2);
            plot([0 L(1) L(1) 0 0], [0 0 L(2) L(2) 0],'k-','linewidth',2);
        %}
    end

    % Step 2. adjust simulation based on net motion of system
    if PFIRE > 0
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
        V = V - 0.5*dt*Pi;

        % reset velocities to 0
        vx = zeros(n,2);
        vy = zeros(n,2);
        Pi = 0.0;
    end

    % update vnorm and fnorm for FIRE
    vnorm = sqrt(sum(vx.*vx,'all') + sum(vy.*vy,'all') + Pi*Pi);
    fnorm = sqrt(sum((fx - 0.5*vx*(Pi/V)).^2,'all') + sum((fy - 0.5*vy*(Pi/V)).^2,'all') + (P-P0)^2);

    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*((fx - 0.5*vx*(Pi/V))./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*((fy - 0.5*vy*(Pi/V))./fnorm)*vnorm;
        Pi = (1 - alpha)*Pi + alpha*((P - P0)/fnorm)*vnorm;
    end

    % do first verlet update for vertices (assume unit mass)
    x = x + dt*(vx + 0.5*x*(Pi/V));
    y = y + dt*(vy + 0.5*y*(Pi/V));
    V = V + dt*Pi;
    Lbase = sqrt(V/boxAR);
    Lxold = L(1);
    Lyold = L(2);
    L(1) = Lbase;
    L(2) = boxAR*Lbase;
    
    % scale position of box edge pins
    for nn = 1:2
        for ii = 1:n
            bwtmp = bw(ii,1,nn);
            if bwtmp > 0
                bw(ii,2,nn) = (L(1)/Lxold)*bw(ii,2,nn);
                bw(ii,3,nn) = (L(2)/Lyold)*bw(ii,3,nn);
            end
        end
    end
    
    % figure out forces on all vertices for both cells
    fx = zeros(n,2);
    fy = zeros(n,2);
    for nn = 1:2
        % shape parameters
        xtmp = x(:,nn);
        ytmp = y(:,nn);
        rtmp = r(:,nn);
        a0tmp = a0(nn);
        l0tmp = l0(:,nn);
        th0tmp = th0(:,nn);

        % -- perimeter force

        % segment vectors
        lvx = xtmp(ip1) - xtmp;
        lvy = ytmp(ip1) - ytmp;

        % update perimeter segment lengths
        l = sqrt(lvx.^2 + lvy.^2);

        % segment unit vectos
        ulvx = lvx./l;
        ulvy = lvy./l;

        % segment strain
        dli = ((l./l0tmp) - 1.0)./l0tmp;
        dlim1 = ((l(im1)./l0tmp(im1)) - 1.0)./l0tmp(im1);

        % perimeter force at this iteration
        flx = kl*(dli.*ulvx - dlim1.*ulvx(im1));
        fly = kl*(dli.*ulvy - dlim1.*ulvy(im1));

        % add to total force, potential
        fx(:,nn) = fx(:,nn) + flx;
        fy(:,nn) = fy(:,nn) + fly;

        % -- area force
        % update area
        atmp = polyarea(xtmp, ytmp);
        da = (atmp/a0tmp) - 1.0;

        % area force at this iteration
        fax = 0.5*da.*(ytmp(im1) - ytmp(ip1))./a0tmp;
        fay = 0.5*da.*(xtmp(ip1) - xtmp(im1))./a0tmp;

        % add to total force, potential
        fx(:,nn) = fx(:,nn) + fax;
        fy(:,nn) = fy(:,nn) + fay;
        
        % -- bending force

        % get sine + cosine
        si = lvx.*lvy(im1) - lvy.*lvx(im1);
        ci = lvx.*lvx(im1) + lvy.*lvy(im1);
        thi = atan2(si,ci);
        dth = thi - th0tmp;

        % get normal vector
        nix = lvy;
        niy = -lvx;

        % construct force components
        fbix = (dth - dth(ip1)).*nix./(l.^2);
        fbiy = (dth - dth(ip1)).*niy./(l.^2);

        fbim1x = -fbix(im1);
        fbim1y = -fbiy(im1);

        % add up forces
        fbx = kb*(fbim1x + fbix);
        fby = kb*(fbim1y + fbiy);

        % add to force, potential
        fx(:,nn) = fx(:,nn) + fbx;
        fy(:,nn) = fy(:,nn) + fby;
        
        
        % check wall forces
        fx_wall = zeros(n,1);
        fy_wall = zeros(n,1);
        
        % wall indexes
%         left_idx = xtmp < rtmp;
%         bottom_idx = ytmp < rtmp;
%         right_idx = xtmp > L(1) - rtmp;
%         top_idx = ytmp > L(2) - rtmp;
        left_idx = xtmp < rtmp | bw(:,1,nn) == 1;
        bottom_idx = ytmp < rtmp | bw(:,1,nn) == 2;
        right_idx = xtmp > L(1) - rtmp | bw(:,1,nn) == 3;
        top_idx = ytmp > L(2) - rtmp | bw(:,1,nn) == 4;
                
        % get wall forces
        fx_wall(left_idx) = (kw./rtmp(left_idx)).*(1.0 - (xtmp(left_idx)./rtmp(left_idx)));
        fy_wall(bottom_idx) = (kw./rtmp(bottom_idx)).*(1.0 - (ytmp(bottom_idx)./rtmp(bottom_idx)));
        
        fx_wall(right_idx) = -(kw./rtmp(right_idx)).*(1.0 - ((L(1) - xtmp(right_idx))./rtmp(right_idx)));
        fy_wall(top_idx) = -(kw./rtmp(top_idx)).*(1.0 - ((L(2) - ytmp(top_idx))./rtmp(top_idx)));
        
        fx(:,nn) = fx(:,nn) + fx_wall;
        fy(:,nn) = fy(:,nn) + fy_wall;
        
        % update bw labels
        ov_left = xtmp < rtmp & bw(:,1,nn) == 0;
        ov_bottom = ytmp < rtmp & bw(:,1,nn) == 0;
        ov_right = xtmp > L(1) - rtmp & bw(:,1,nn) == 0;
        ov_top = ytmp > L(2) - rtmp & bw(:,1,nn) == 0;
        
        bw(ov_left,1,nn) = -1*ones(sum(ov_left),1);
        bw(ov_bottom,1,nn) = -1*ones(sum(ov_bottom),1);
        bw(ov_right,1,nn) = -1*ones(sum(ov_right),1);
        bw(ov_top,1,nn) = -1*ones(sum(ov_top),1);
        
        
        % check contact forces
        if nn == 1
            for ii = 1:n
                xi = xtmp(ii);
                yi = ytmp(ii);
                ri = rtmp(ii);
                for jj = 1:n
                    % contact distance
                    sij = ri + r(jj,2);
                    
                    % distance
                    xj = x(jj,2);
                    yj = y(jj,2);
                    
                    dx = xj - xi;
                    dy = yj - yi;
                    
                    rij = sqrt(dx*dx + dy*dy);
                    
                    % check contacts
                    if (rij < sij || bij(ii,jj) == 1)
                        ftmp = (kc/sij)*(1 - (rij/sij));
                        fxtmp = ftmp*(dx/rij);
                        fytmp = ftmp*(dy/rij);
                        
                        fx(ii,1) = fx(ii,1) - fxtmp;
                        fy(ii,1) = fy(ii,1) - fytmp;
                        
                        fx(jj,2) = fx(jj,2) + fxtmp;
                        fy(jj,2) = fy(jj,2) + fytmp;
                    end
                end
                
%                 % also check wall bonds
%                 if bw(ii,1,nn) > 0
%                     dx = xi - bw(ii,2,nn);
%                     dy = yi - bw(ii,3,nn);
%                     dr = sqrt(dx*dx + dy*dy);
%                     ftmp = (kw/ri)*(1 - (dr/ri));
%                     fxwb = ftmp*(dx/dr);
%                     fywb = ftmp*(dy/dr);
%                     fx(ii,nn) = fx(ii,nn) + fxwb;
%                     fy(ii,nn) = fy(ii,nn) + fywb;
%                 end
            end
        end
    end
    
    % update pressure based on positions
    Pext = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
    Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
    P = Pext + Pint;

    % step 3 (VV)
    vx = vx + 0.5*dt*(fx - 0.5*vx*(Pi/V));
    vy = vy + 0.5*dt*(fy - 0.5*vy*(Pi/V));
    Pi = Pi + 0.5*dt*(P - P0)/PM;

    % update force check
    fcheck = sqrt(sum(sum(fx(:).^2) + sum(fy(:).^2) + (P-P0)^2)/(2.0*n));
end



end


function P = mesoDimerPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw)
%% FUNCTION to measure pressure of meso dimer

% number of vertices
n = size(x,1);

% indexing
ip1 = [2:n 1];

% loop
Lbase = L(1);
alpha = L(2)/L(1);
P = 0.0;
for nn = 1:2
    % shape parameters
    xtmp = x(:,nn);
    ytmp = y(:,nn);
    rtmp = r(:,nn);
    a0tmp = a0(nn);
    l0tmp = l0(:,nn);
    
    % -- perimeter contribution

    % segment vectors
    lvx = xtmp(ip1) - xtmp;
    lvy = ytmp(ip1) - ytmp;

    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    dli = (l./l0tmp) - 1.0;
    
    % add to pressure
    P = P + (kl/Lbase)*sum((l./l0tmp).*dli);

    % -- area contribution
    
    % polygon area
    atmp = polyarea(xtmp, ytmp);
    da = (atmp/a0tmp) - 1.0;
    
    % add to pressure
    P = P + (2.0/Lbase)*((atmp/a0tmp).*da);
    
    % -- wall contributions
    dUwdL = 0.0;
    for vv = 1:n
        xv = xtmp(vv);
        yv = ytmp(vv);
        rv = rtmp(vv);
        
        % x walls
        if xv < rv || bw(vv,1,nn) == 1
            ftmp = kw*(1 - (xv/rv))/rv;
            dUwdL = dUwdL - ftmp*(xv/Lbase);
        elseif xv > L(1) - rv || bw(vv,1,nn) == 3
            ftmp = kw*(1 - ((L(1) - xv)/rv))/rv;
            dUwdL = dUwdL + ftmp*((xv/Lbase) - 1);
        end
        
        % y walls
        if yv < rv || bw(vv,1,nn) == 2
            ftmp = kw*(1 - (yv/rv))/rv;
            dUwdL = dUwdL - ftmp*(yv/Lbase);
        elseif yv > L(2) - rv || bw(vv,1,nn) == 4
            ftmp = kw*(1 - ((L(2) - yv)/rv))/rv;
            dUwdL = dUwdL + ftmp*((yv/Lbase) - alpha);
        end
        
%         % check bonds
%         if bw(vv,1,nn) > 0
%             dx = xv - bw(vv,2,nn);
%             dy = yv - bw(vv,3,nn);
%             dr = sqrt(dx*dx + dy*dy);
%             ftmp = -(kw/rv)*(1 - (dr/rv));
%             dUwdL = dUwdL + ftmp*(dr/Lbase);
%         end
    end
    P = P + dUwdL;
    
    
    % -- interaction contributions
    if nn == 1
        for ii = 1:n
            xi = xtmp(ii);
            yi = ytmp(ii);
            for jj = 1:n
                % contact distance
                sij = rtmp(ii) + r(jj,2);

                % distance
                xj = x(jj,2);
                yj = y(jj,2);

                dx = xj - xi;
                dy = yj - yi;

                rij = sqrt(dx*dx + dy*dy);

                % check contacts
                if (rij < sij || bij(ii,jj) == 1)
                    ftmp = (kc/sij)*(1 - (rij/sij));
                    P = P - ftmp*(rij/Lbase);
                end
            end
        end
    end
end

% scale 
P = -P/(2.0*alpha*Lbase);




end
