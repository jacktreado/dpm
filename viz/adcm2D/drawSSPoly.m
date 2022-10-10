function drawSSPoly(x, y, r, clr, Lx, Ly, pbcs)
% number of vertices
n = length(x);

% segments
ip1 = [2:n 1];
lx = x(ip1) - x;
ly = y(ip1) - y;
l = sqrt(lx.^2 + ly.^2);
th = atan2(ly, lx) * (180 / pi);

if pbcs == 1
    dbox = -1:1;
else
    dbox = 0;
end

% draw vertices as circles
for ii = 1:n
    for xx = dbox
        for yy = dbox
            rectangle('Position', [x(ii) - r(ii) + xx*Lx, y(ii) - r(ii) + Ly*yy, 2 * r(ii), 2 * r(ii)], 'Curvature', [1 1], 'FaceColor', 'w', 'EdgeColor', clr, 'LineWidth', 2);
        end
    end
end

% draw patch
for xx = dbox
    for yy = dbox
        patch('Faces',[1:n 1],'Vertices',[x + xx*Lx y + yy*Ly],'FaceColor',clr,'EdgeColor','k');
    end
end

% loop over segments
for ii = 1:n
    % plot flat segment starting at x(ii) with length l(ii)
    % NOTE: unrotated segment is constructed with y - r, not y + r, because
    % rotation must be on top half of segment (it's confusing)
    for xx = dbox
        for yy = dbox
            p = patch([x(ii); x(ii) + l(ii); x(ii) + l(ii); x(ii)] + Lx*xx,[y(ii); y(ii); y(ii) - r(ii); y(ii) - r(ii)] + Ly*yy, 'k', 'EdgeColor', clr, 'FaceColor', 'w', 'LineWidth', 2);
            rotate(p,[0 0 1],th(ii),[x(ii) + Lx*xx ,y(ii) + Ly*yy,0]);
        end
    end
    
end
end