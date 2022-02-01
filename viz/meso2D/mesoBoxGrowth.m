function mesoBoxGrowth(savestr,N,NGROWTH,calA0_base,kl_base,kb_base,kc,dl0,da0,cB,P0)
%% FUNCTION to run mesoBoxGrowth in single, self-contained function, save output

% constants
phi0        = 0.1;
Pcmp        = 1e-2;
bbreak      = 1.1;
n           = 32;
l0_base     = sqrt(4.0*pi*calA0_base)/n;
kw          = kc;


% shape parameters
a0          = ones(1,N);
l0          = l0_base*ones(n,N);
th0         = zeros(n,N);
kl          = kl_base*ones(n,N);
kb          = kb_base*ones(n,N);
kbinit      = zeros(n,N);
r           = zeros(n,N);

% arrays for degrees of freedom and vertex sizes
x           = zeros(n,N);
y           = zeros(n,N);

% bonded matrices
NVTOT       = N*n;
bij         = zeros(NVTOT,NVTOT);
bw          = zeros(n,3,N);


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

% initial box length (as square, can be edited)
L0 = sqrt(N/phi0);
L = L0*ones(2,1);

% set positions of vertices (around center of box)
for nn = 1:N
    x(:,nn) = 0.5*L0 + 2*r0*cos((2.0*pi*(nn-1))/N) + x0;
    y(:,nn) = 0.5*L0 + 2*r0*sin((2.0*pi*(nn-1))/N) + y0;
    r(:,nn) = repmat(0.5*l0_base,n,1);
end


% do pressure minimization
fprintf('** Initial H min. at Pcmp=%0.5g\n',Pcmp);
[x,y,fx,fy,L,bw,bij] = mesoBoxHMin(x,y,r,a0,l0,th0,Pcmp,L,kl,kbinit,kc,kw,bij,bw);

Pext = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
P = Pext + Pint;
fprintf('** After H min., pressure P = %0.5g, Pcmp = %0.5g\n',P,Pcmp);



%% Set up bond network, relax to production pressure

% initialize pairwise bonds
for nn = 1:N
    for ii = 1:n
        xi = x(ii,nn);
        yi = y(ii,nn);
        ri = r(ii,nn);
        for mm = (nn+1):N
            for jj = 1:n
                xj = x(jj,mm);
                yj = y(jj,mm);
                rj = r(jj,mm);
                sij = ri + rj;
                dx = xj - xi;
                dy = yj - yi;
                rij = sqrt(dx*dx + dy*dy);
                
                % global vertex indices
                gi = (nn-1)*n + ii;
                gj = (mm-1)*n + jj;
                if rij < sij
                    bij(gi,gj) = 1;
                end
            end
        end
    end
end

% initialize wall bonds
for nn = 1:N
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


% do pressure minimization with bonds
fprintf('** Initial H min. WITH BONDS at P0=%0.5g\n',P0);
[x,y,fx,fy,L,bw,bij] = mesoBoxHMin(x,y,r,a0,l0,th0,P0,L,kl,kbinit,kc,kw,bij,bw);

Pext = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
P = Pext + Pint;
fprintf('** After H min. with bonds, pressure P = %0.5g, P0 = %0.5g\n',P,P0);

% set preferred angles to current setup
thi = currentAngles(x,y);
th0 = thi;


%% Grow cells

% save during growth steps
xList = cell(NGROWTH,1);
yList = cell(NGROWTH,1);
rList = cell(NGROWTH,1);
LList = zeros(NGROWTH,2);
aList = zeros(NGROWTH,N);
a0List = zeros(NGROWTH,N);
lList = cell(NGROWTH,1);
l0List = cell(NGROWTH,1);
thiList = cell(NGROWTH,1);
th0List = cell(NGROWTH,1);
calAList = zeros(NGROWTH,N);
calA0List = zeros(NGROWTH,N);
bijList = cell(NGROWTH,1);
bwList = cell(NGROWTH,1);


for gg = 1:NGROWTH
    fprintf('Enthalpy min + growth, step gg = %d\n',gg);

    % grow area, void perimeter
    for nn = 1:N
        a0(nn) = a0(nn)*(1 + dl0*da0)^2;
        r(:,nn) = r(:,nn)*(1 + dl0*da0);
        for ii = 1:n
            gi = (nn-1)*n + ii;
            im1 = mod(ii+n-2,n)+1;
            if bw(ii,1,nn) == 0 && sum(bij(gi,:)) == 0
                l0([im1 ii],nn) = l0([im1 ii],nn)*(1 + dl0);
                th0([im1 ii],nn) = th0([im1 ii],nn) + dl0*cB;
            elseif sum(bij(gi,:)) ~= 0 || bw(ii,1,nn) ~= 0
                th0(ii,nn) = th0(ii,nn)*(1.0 - cB*dl0);
            end
        end
    end
                
    % enthalpy minimize
    fprintf('Minimizing H with void perimeter growth\n');
    [x,y,fx,fy,L,bw,bij] = mesoBoxHMin(x,y,r,a0,l0,th0,P0,L,kl,kb,kc,kw,bij,bw);
    
    % compute pressure, print
    Pext = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
    Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
    P = Pext + Pint;
    lx = x([2:n 1],:) - x;
    ly = y([2:n 1],:) - y;
    l = sqrt(lx.^2 + ly.^2);
    a = polyarea(x,y);
    da = (a./a0) - 1;
    p = sum(l,1);
    calA = p.^2./(4.0*pi*a);
    calA0 = sum(l0).^2./(4.0*pi*a0);
    thi = currentAngles(x,y);
    
    fprintf('\n** Enthalpy minimized + printing\n');
    fprintf('** Pext = %0.5g\n',Pext);
    fprintf('** Pint = %0.5g\n',Pint);
    fprintf('** P = %0.5g\n',P);
    fprintf('** L = %0.5g\n',L(1));
    fprintf('** calA1 = %0.5g, calA2 = %0.5g\n',calA(1),calA(2));
    fprintf('** calA01 = %0.5g, calA02 = %0.5g\n',calA0(1),calA0(2));
    fprintf('** phi = %0.5g\n',sum(a)/(L(1)*L(2)));
    fprintf('** area strain = %0.5g, %0.5g\n\n\n',da(1),da(2));
    
    % save configuration
    xList{gg} = x;
    yList{gg} = y;
    rList{gg} = r;
    LList(gg,:) = L;
    aList(gg,:) = a;
    a0List(gg,:) = a0;
    lList{gg} = l;
    l0List{gg} = l0;
    thiList{gg} = thi;
    th0List{gg} = th0;
    calAList(gg,:) = calA;
    calA0List(gg,:) = calA0;
    bijList{gg} = bij;
    bwList{gg} = bw;
    
    
    % break long bonds
    for nn = 1:N
        for ii = 1:n
            xi = x(ii,nn);
            yi = y(ii,nn);
            ri = r(ii,nn);
            gi = (nn-1)*n + ii;
            for mm = (nn+1):N
                for jj = 1:n
                    gj = (mm-1)*n + jj;
                    if bij(gi,gj) == 1
                        xj = x(jj,mm);
                        yj = y(jj,mm);
                        rj = r(jj,mm);
                        sij = ri + rj;
                        dx = xj - xi;
                        dy = yj - yi;
                        rij = sqrt(dx*dx + dy*dy);
                        if rij > bbreak*sij
                           bij(gi,gj) = 0;
                           bij(gj,gi) = 0; 
                        end
                    end
                end
            end
        end
    end

    % break wall bonds
    for nn = 1:N
        for ii = 1:n
            % wall bond
            bwtmp = bw(ii,1,nn);
            if bwtmp > 0
                dx = x(ii,nn) - bw(ii,2,nn);
                dy = y(ii,nn) - bw(ii,3,nn);
                dr = sqrt(dx*dx + dy*dy);
                if dr > bbreak*r(ii,nn)
                    bw(ii,1,nn) = 0;
                end
            end
        end
    end
end

% save
save(savestr,...
    'n','N','NGROWTH','calA0_base','kl_base',...
    'kb_base','kc','dl0','da0','cB','P0',...
    'xList','yList','rList','LList',...
    'aList','a0List','lList','l0List',...
    'thiList','th0List','calAList','calA0List',...
    'bijList','bwList');

end




%% Enthalpy minimization function

function [x,y,fx,fy,L,bw,bij] = mesoBoxHMin(x,y,r,a0,l0,th0,P0,L,kl,kb,kc,kw,bij,bw)
%% FUNCTION to relax total potential energy using fire for meso dimer

% number of cells and vertices
[n, N] = size(x);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];

% forces and velocities
vx = zeros(n,N);
vy = zeros(n,N);
fx = zeros(n,N);
fy = zeros(n,N);

% enthalpy variables
V = L(1)*L(2);
Pi = 0.0;
P = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
boxAR = L(2)/L(1);
PM = 1;

% FIRE parameters
alpha0      = 0.2;
finc        = 1.1;
fdec        = 0.5;
falpha      = 0.99;
dt0         = 0.005;
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
    if mod(it,5000) == 0
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
        % &&&&&&&&&&&&&&&&&&&&&&
        % Plotting for debugging
        % &&&&&&&&&&&&&&&&&&&&&&

        figure(1), clf, hold on, box on;
        for nn = 1:N
            for vv = 1:n
                rectangle('Position',[x(vv,nn)-r(vv,nn),y(vv,nn)-r(vv,nn),2.0*r(vv,nn),2.0*r(vv,nn)],'Curvature',[1 1],'FaceColor','b');
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
        vx = zeros(n,N);
        vy = zeros(n,N);
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
    for nn = 1:N
        for ii = 1:n
            bwtmp = bw(ii,1,nn);
            if bwtmp > 0
                bw(ii,2,nn) = (L(1)/Lxold)*bw(ii,2,nn);
                bw(ii,3,nn) = (L(2)/Lyold)*bw(ii,3,nn);
            end
        end
    end
    
    % figure out forces on all vertices for both cells
    fx = zeros(n,N);
    fy = zeros(n,N);
    for nn = 1:N
        % shape parameters
        xtmp = x(:,nn);
        ytmp = y(:,nn);
        rtmp = r(:,nn);
        a0tmp = a0(nn);
        l0tmp = l0(:,nn);
        th0tmp = th0(:,nn);
        kltmp = kl(:,nn);
        kbtmp = kb(:,nn);

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
        flx = kltmp.*dli.*ulvx - kltmp(im1).*dlim1.*ulvx(im1);
        fly = kltmp.*dli.*ulvy - kltmp(im1).*dlim1.*ulvy(im1);

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
        fbx = kbtmp(im1).*fbim1x + kbtmp.*fbix;
        fby = kbtmp(im1).*fbim1y + kbtmp.*fbiy;

        % add to force, potential
        fx(:,nn) = fx(:,nn) + fbx;
        fy(:,nn) = fy(:,nn) + fby;
        
        
        % check wall forces
        fx_wall = zeros(n,1);
        fy_wall = zeros(n,1);
        
        % wall indexes
        left_idx = xtmp < rtmp;
        bottom_idx = ytmp < rtmp;
        right_idx = xtmp > L(1) - rtmp;
        top_idx = ytmp > L(2) - rtmp;
%         left_idx = xtmp < rtmp | bw(:,1,nn) == 1;
%         bottom_idx = ytmp < rtmp | bw(:,1,nn) == 2;
%         right_idx = xtmp > L(1) - rtmp | bw(:,1,nn) == 3;
%         top_idx = ytmp > L(2) - rtmp | bw(:,1,nn) == 4;
                
        % get wall forces
        fx_wall(left_idx) = (kw./rtmp(left_idx)).*(1.0 - (xtmp(left_idx)./rtmp(left_idx)));
        fy_wall(bottom_idx) = (kw./rtmp(bottom_idx)).*(1.0 - (ytmp(bottom_idx)./rtmp(bottom_idx)));
        
        fx_wall(right_idx) = -(kw./rtmp(right_idx)).*(1.0 - ((L(1) - xtmp(right_idx))./rtmp(right_idx)));
        fy_wall(top_idx) = -(kw./rtmp(top_idx)).*(1.0 - ((L(2) - ytmp(top_idx))./rtmp(top_idx)));
        
        fx(:,nn) = fx(:,nn) + fx_wall;
        fy(:,nn) = fy(:,nn) + fy_wall;
        
        % update wall bonds for growth
        ov_left = xtmp < rtmp & bw(:,1,nn) == 0;
        ov_bottom = ytmp < rtmp & bw(:,1,nn) == 0;
        ov_right = xtmp > L(1) - rtmp & bw(:,1,nn) == 0;
        ov_top = ytmp > L(2) - rtmp & bw(:,1,nn) == 0;
        
        off_left = xtmp > rtmp & bw(:,1,nn) == -1;
        off_bottom = ytmp > rtmp & bw(:,1,nn) == -1;
        off_right = xtmp < L(1) - rtmp & bw(:,1,nn) == -1;
        off_top = ytmp < L(2) - rtmp & bw(:,1,nn) == -1;

        bw(ov_left,1,nn) = -1*ones(sum(ov_left),1);
        bw(ov_bottom,1,nn) = -1*ones(sum(ov_bottom),1);
        bw(ov_right,1,nn) = -1*ones(sum(ov_right),1);
        bw(ov_top,1,nn) = -1*ones(sum(ov_top),1);
        
        bw(off_left,1,nn) = zeros(sum(off_left),1);
        bw(off_bottom,1,nn) = zeros(sum(off_bottom),1);
        bw(off_right,1,nn) = zeros(sum(off_right),1);
        bw(off_top,1,nn) = zeros(sum(off_top),1);
        
        % -- interaction contributions
        for ii = 1:n
            xi = xtmp(ii);
            yi = ytmp(ii);
            ri = rtmp(ii);
            for mm = (nn+1):N
                for jj = 1:n
                    % contact distance
                    sij = ri + r(jj,mm);

                    % distance
                    xj = x(jj,mm);
                    yj = y(jj,mm);

                    dx = xj - xi;
                    dy = yj - yi;

                    rij = sqrt(dx*dx + dy*dy);
                    
                    % global vertex indices
                    gi = (nn-1)*n + ii;
                    gj = (mm-1)*n + jj;
                    
                    % check contacts
                    if (rij < sij || bij(gi,gj) == 1)
                        ftmp = (kc/sij)*(1 - (rij/sij));
                        fxtmp = ftmp*(dx/rij);
                        fytmp = ftmp*(dy/rij);
                        
                        fx(ii,nn) = fx(ii,nn) - fxtmp;
                        fy(ii,nn) = fy(ii,nn) - fytmp;
                        
                        fx(jj,mm) = fx(jj,mm) + fxtmp;
                        fy(jj,mm) = fy(jj,mm) + fytmp;
                        
                        if bij(gi,gj) == 0
                            bij(gi,gj) = -1;
                            bij(gj,gi) = bij(gi,gj);
                        end
                    else
                        if bij(gi,gj) == -1
                            bij(gi,gj) = 0;
                            bij(gj,gi) = bij(gi,gj);
                        end
                    end
                end
            end

            % also check wall bonds
            if bw(ii,1,nn) > 0
                dx = xi - bw(ii,2,nn);
                dy = yi - bw(ii,3,nn);
                dr = sqrt(dx*dx + dy*dy);
                ftmp = (kw/ri)*(1 - (dr/ri));
                fxwb = ftmp*(dx/dr);
                fywb = ftmp*(dy/dr);
                fx(ii,nn) = fx(ii,nn) + fxwb;
                fy(ii,nn) = fy(ii,nn) + fywb;
            end
        end
    end
    
    % update pressure based on positions
    Pext = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw);
    Pint = (sum(fx.*x,'all') + sum(fy.*y,'all'))./(2.0*L(1)*L(2));
    P = Pext + Pint;

    % step 3 (VV)
    vx = vx + 0.5*dt*(fx - 0.5*vx*(Pi/V));
    vy = vy + 0.5*dt*(fy - 0.5*vy*(Pi/V));
    Pi = Pi + 0.5*dt*(P - P0)/PM;

    % update force check
    fcheck = sqrt(sum(sum(fx(:).^2) + sum(fy(:).^2) + (P-P0)^2)/(N*n));
end


end









%% Energy minimization function

function [x,y,fx,fy] = mesoBoxEMin(x,y,r,a0,l0,th0,L,kl,kb,kc,kw,bij,bw)
%% FUNCTION to relax total potential energy using fire for meso dimer

% number of cells and vertices
[n, N] = size(x);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];

% forces and velocities
vx = zeros(n,N);
vy = zeros(n,N);
fx = zeros(n,N);
fy = zeros(n,N);

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
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % Step 1. calculate P, fnorm, vnorm
    PFIRE = sum(vx.*fx,'all') + sum(vy.*fy,'all');

    % plot FIRE information 
    if mod(it,5000) == 0
        fprintf('\n\t && On FIRE EMin step %d\n',it);
        fprintf('\t ** F = %0.5g\n',fcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** PFIRE = %0.5g\n',PFIRE);
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** npPMin = %d\n',npPMin);

        
        %{
        % &&&&&&&&&&&&&&&&&&&&&&
        % Plotting for debugging
        % &&&&&&&&&&&&&&&&&&&&&&

        figure(1), clf, hold on, box on;
        for nn = 1:N
            for vv = 1:n
                rectangle('Position',[x(vv,nn)-r(vv,nn),y(vv,nn)-r(vv,nn),2.0*r(vv,nn),2.0*r(vv,nn)],'Curvature',[1 1],'FaceColor','b');
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

        % reset velocities to 0
        vx = zeros(n,N);
        vy = zeros(n,N);
    end

    % update vnorm and fnorm for FIRE
    vnorm = sqrt(sum(vx.*vx,'all') + sum(vy.*vy,'all'));
    fnorm = sqrt(sum(fx,'all') + sum(fy,'all'));

    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*(fx./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*(fy./fnorm)*vnorm;
    end

    % do first verlet update for vertices (assume unit mass)
    x = x + dt*vx;
    y = y + dt*vy;
    
    % figure out forces on all vertices for both cells
    fx = zeros(n,N);
    fy = zeros(n,N);
    for nn = 1:N
        % shape parameters
        xtmp = x(:,nn);
        ytmp = y(:,nn);
        rtmp = r(:,nn);
        a0tmp = a0(nn);
        l0tmp = l0(:,nn);
        th0tmp = th0(:,nn);
        kltmp = kl(:,nn);
        kbtmp = kb(:,nn);

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
        flx = kltmp.*dli.*ulvx - kltmp(im1).*dlim1.*ulvx(im1);
        fly = kltmp.*dli.*ulvy - kltmp(im1).*dlim1.*ulvy(im1);

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
        fbx = kbtmp(im1).*fbim1x + kbtmp.*fbix;
        fby = kbtmp(im1).*fbim1y + kbtmp.*fbiy;

        % add to force, potential
        fx(:,nn) = fx(:,nn) + fbx;
        fy(:,nn) = fy(:,nn) + fby;
        
        
        % check wall forces
        fx_wall = zeros(n,1);
        fy_wall = zeros(n,1);
        
        % wall indexes
        left_idx = xtmp < rtmp;
        bottom_idx = ytmp < rtmp;
        right_idx = xtmp > L(1) - rtmp;
        top_idx = ytmp > L(2) - rtmp;
%         left_idx = xtmp < rtmp | bw(:,1,nn) == 1;
%         bottom_idx = ytmp < rtmp | bw(:,1,nn) == 2;
%         right_idx = xtmp > L(1) - rtmp | bw(:,1,nn) == 3;
%         top_idx = ytmp > L(2) - rtmp | bw(:,1,nn) == 4;
                
        % get wall forces
        fx_wall(left_idx) = (kw./rtmp(left_idx)).*(1.0 - (xtmp(left_idx)./rtmp(left_idx)));
        fy_wall(bottom_idx) = (kw./rtmp(bottom_idx)).*(1.0 - (ytmp(bottom_idx)./rtmp(bottom_idx)));
        
        fx_wall(right_idx) = -(kw./rtmp(right_idx)).*(1.0 - ((L(1) - xtmp(right_idx))./rtmp(right_idx)));
        fy_wall(top_idx) = -(kw./rtmp(top_idx)).*(1.0 - ((L(2) - ytmp(top_idx))./rtmp(top_idx)));
        
        fx(:,nn) = fx(:,nn) + fx_wall;
        fy(:,nn) = fy(:,nn) + fy_wall;
        
        % -- interaction contributions
        for ii = 1:n
            xi = xtmp(ii);
            yi = ytmp(ii);
            ri = rtmp(ii);
            for mm = (nn+1):N
                for jj = 1:n
                    % contact distance
                    sij = ri + r(jj,mm);

                    % distance
                    xj = x(jj,mm);
                    yj = y(jj,mm);

                    dx = xj - xi;
                    dy = yj - yi;

                    rij = sqrt(dx*dx + dy*dy);
                    
                    % global vertex indices
                    gi = (nn-1)*n + ii;
                    gj = (mm-1)*n + jj;
                    
                    % check contacts
                    if (rij < sij || bij(gi,gj) == 1)
                        ftmp = (kc/sij)*(1 - (rij/sij));
                        fxtmp = ftmp*(dx/rij);
                        fytmp = ftmp*(dy/rij);
                        
                        fx(ii,nn) = fx(ii,nn) - fxtmp;
                        fy(ii,nn) = fy(ii,nn) - fytmp;
                        
                        fx(jj,mm) = fx(jj,mm) + fxtmp;
                        fy(jj,mm) = fy(jj,mm) + fytmp;
                    end
                end
            end
            
            % also check wall bonds
            if bw(ii,1,nn) > 0
                dx = xi - bw(ii,2,nn);
                dy = yi - bw(ii,3,nn);
                dr = sqrt(dx*dx + dy*dy);
                ftmp = (kw/ri)*(1 - (dr/ri));
                fxwb = ftmp*(dx/dr);
                fywb = ftmp*(dy/dr);
                fx(ii,nn) = fx(ii,nn) + fxwb;
                fy(ii,nn) = fy(ii,nn) + fywb;
            end
        end
    end

    % step 3 (VV)
    vx = vx + 0.5*dt*fx;
    vy = vy + 0.5*dt*fy;

    % update force check
    fcheck = sqrt(sum(sum(fx(:).^2) + sum(fy(:).^2))/(N*n));
end


end










%% External contribution to pressure

function P = mesoBoxPressure(x,y,r,a0,l0,L,kl,kc,kw,bij,bw)
%% FUNCTION to measure pressure of meso dimer

% number of cells and vertices
[n, N] = size(x);

% indexing
ip1 = [2:n 1];

% loop
Lbase = L(1);
alpha = L(2)/L(1);
P = 0.0;
for nn = 1:N
    % shape parameters
    xtmp = x(:,nn);
    ytmp = y(:,nn);
    rtmp = r(:,nn);
    a0tmp = a0(nn);
    l0tmp = l0(:,nn);
    kltmp = kl(:,nn);
    
    % -- perimeter contribution

    % segment vectors
    lvx = xtmp(ip1) - xtmp;
    lvy = ytmp(ip1) - ytmp;

    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    dli = (l./l0tmp) - 1.0;
    
    % add to pressure
    P = P + (1.0/Lbase)*sum(((kltmp.*l)./l0tmp).*dli);

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
        
%         % -- no wall bonds
%         % x walls
%         if xv < rv || bw(vv,1,nn) == 1
%             ftmp = kw*(1 - (xv/rv))/rv;
%             dUwdL = dUwdL - ftmp*(xv/Lbase);
%         elseif xv > L(1) - rv || bw(vv,1,nn) == 3
%             ftmp = kw*(1 - ((L(1) - xv)/rv))/rv;
%             dUwdL = dUwdL + ftmp*((xv/Lbase) - 1);
%         end
%         
%         % y walls
%         if yv < rv || bw(vv,1,nn) == 2
%             ftmp = kw*(1 - (yv/rv))/rv;
%             dUwdL = dUwdL - ftmp*(yv/Lbase);
%         elseif yv > L(2) - rv || bw(vv,1,nn) == 4
%             ftmp = kw*(1 - ((L(2) - yv)/rv))/rv;
%             dUwdL = dUwdL + ftmp*((yv/Lbase) - alpha);
%         end
        
        % -- yes wall bonds
        % x walls
        if xv < rv
            ftmp = kw*(1 - (xv/rv))/rv;
            dUwdL = dUwdL - ftmp*(xv/Lbase);
        elseif xv > L(1) - rv
            ftmp = kw*(1 - ((L(1) - xv)/rv))/rv;
            dUwdL = dUwdL + ftmp*((xv/Lbase) - 1);
        end
        
        % y walls
        if yv < rv
            ftmp = kw*(1 - (yv/rv))/rv;
            dUwdL = dUwdL - ftmp*(yv/Lbase);
        elseif yv > L(2) - rv
            ftmp = kw*(1 - ((L(2) - yv)/rv))/rv;
            dUwdL = dUwdL + ftmp*((yv/Lbase) - alpha);
        end
        
        % check bonds
        if bw(vv,1,nn) > 0
            dx = xv - bw(vv,2,nn);
            dy = yv - bw(vv,3,nn);
            dr = sqrt(dx*dx + dy*dy);
            ftmp = -(kw/rv)*(1 - (dr/rv));
            dUwdL = dUwdL + ftmp*(dr/Lbase);
        end
    end
    P = P + dUwdL;
    
    
    % -- interaction contributions
    for ii = 1:n
        xi = xtmp(ii);
        yi = ytmp(ii);
        ri = rtmp(ii);
        for mm = (nn+1):N
            for jj = 1:n
                % contact distance
                sij = ri + r(jj,mm);

                % distance
                xj = x(jj,mm);
                yj = y(jj,mm);

                dx = xj - xi;
                dy = yj - yi;

                rij = sqrt(dx*dx + dy*dy);

                % global vertex indices
                gi = (nn-1)*n + ii;
                gj = (mm-1)*n + jj;

                % check contacts
                if (rij < sij || bij(gi,gj) == 1)
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










%% Total potential energy

function U = mesoBoxU(x,y,r,a0,l0,th0,L,kl,kb,kc,kw,bij,bw)
%% FUNCTION to measure pressure of meso dimer

% number of cells and vertices
[n, N] = size(x);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];
U = 0.0;
for nn = 1:N
    % shape parameters
    xtmp = x(:,nn);
    ytmp = y(:,nn);
    rtmp = r(:,nn);
    a0tmp = a0(nn);
    l0tmp = l0(:,nn);
    th0tmp = th0(:,nn);
    kltmp = kl(:,nn);
    kbtmp = kb(:,nn);
    
    % -- perimeter contribution

    % segment vectors
    lvx = xtmp(ip1) - xtmp;
    lvy = ytmp(ip1) - ytmp;

    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    dli = (l./l0tmp) - 1.0;
    
    % add to potential
    U = U + 0.5*sum(kltmp.*dli.^2);

    % -- area contribution
    
    % polygon area
    atmp = polyarea(xtmp, ytmp);
    da = (atmp/a0tmp) - 1.0;
    
    % add to potential
    U = U + 0.5*da^2;
    
    
    % -- bending contribution
    
    % get angles + cosine
    si = lvx.*lvy(im1) - lvy.*lvx(im1);
    ci = lvx.*lvx(im1) + lvy.*lvy(im1);
    thi = atan2(si,ci);
    dth = thi - th0tmp;
    
    % add to potential
    U = U + 0.5*sum(kbtmp.*dth.^2);
    
    % -- wall contributions
    for vv = 1:n
        xv = xtmp(vv);
        yv = ytmp(vv);
        rv = rtmp(vv);
        
        
%         % -- no wall pins
%         % x walls
%         if xv < rv || bw(vv,1,nn) == 1
%             U = U + 0.5*kw*(1 - (xv/rv))^2;
%         elseif xv > L(1) - rv || bw(vv,1,nn) == 3
%             U = U + 0.5*kw*(1 - ((L(1) - xv)/rv))^2;
%         end
%         
%         % y walls
%         if yv < rv || bw(vv,1,nn) == 2
%             U = U + 0.5*kw*(1 - (yv/rv))^2;
%         elseif yv > L(2) - rv || bw(vv,1,nn) == 4
%             U = U + 0.5*kw*(1 - ((L(2) - yv)/rv))^2;
%         end

        % -- with wall pins
        % x walls
        if xv < rv
            U = U + 0.5*kw*(1 - (xv/rv))^2;
        elseif xv > L(1) - rv
            U = U + 0.5*kw*(1 - ((L(1) - xv)/rv))^2;
        end
        
        % y walls
        if yv < rv
            U = U + 0.5*kw*(1 - (yv/rv))^2;
        elseif yv > L(2) - rv
            U = U + 0.5*kw*(1 - ((L(2) - yv)/rv))^2;
        end
        
        % also check wall bonds
        if bw(vv,1,nn) > 0
            dx = xv - bw(vv,2,nn);
            dy = yv - bw(vv,3,nn);
            dr = sqrt(dx*dx + dy*dy);
            U = U + 0.5*kw*(1 - (dr/rv))^2;
        end

    end
    
    
    % -- interaction contributions
    for ii = 1:n
        xi = xtmp(ii);
        yi = ytmp(ii);
        ri = rtmp(ii);
        for mm = (nn+1):N
            for jj = 1:n
                % contact distance
                sij = ri + r(jj,mm);

                % distance
                xj = x(jj,mm);
                yj = y(jj,mm);

                dx = xj - xi;
                dy = yj - yi;

                rij = sqrt(dx*dx + dy*dy);

                % global vertex indices
                gi = (nn-1)*n + ii;
                gj = (mm-1)*n + jj;

                % check contacts
                if (rij < sij || bij(gi,gj) == 1)
                    U = U + 0.5*kc*(1 - (rij/sij))^2;
                end
            end
        end
    end
end
end










%% Return all bending angles

function thi = currentAngles(x,y)
%% FUNCTION to set preferred angles to current angles

% number of cells and vertices
[n, N] = size(x);

% indexing
ip1 = [2:n 1];
im1 = [n 1:n-1];

% loop over cells and vertices
thi = zeros(n,N);
for nn = 1:N
    % shape parameters
    xtmp = x(:,nn);
    ytmp = y(:,nn);
    
    % segment vectors
    lvx = xtmp(ip1) - xtmp;
    lvy = ytmp(ip1) - ytmp;

    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    
    % get sine + cosine
    si = lvx.*lvy(im1) - lvy.*lvx(im1);
    ci = lvx.*lvx(im1) + lvy.*lvy(im1);
    thi = atan2(si,ci);
    
    % assign
    thi(:,nn) = thi;
end
end








