function drawSSPoly(x, y, r, clr)
% number of vertices
n = length(x);

% segments
ip1 = [2:n 1];
lx = x(ip1) - x;
ly = y(ip1) - y;
l = sqrt(lx.^2 + ly.^2);
th = atan2(ly, lx) * (180 / pi);

% draw vertices as circles
for ii = 1:n
    rectangle('Position', [x(ii) - r(ii), y(ii) - r(ii), 2 * r(ii), 2 * r(ii)], 'Curvature', [1 1], 'FaceColor', 'w', 'EdgeColor', clr, 'LineWidth', 2);
end

% draw patch
patch('Faces',[1:n 1],'Vertices',[x y],'FaceColor',clr,'EdgeColor','k');

% loop over segments
for ii = 1:n
    % plot flat segment starting at x(ii) with length l(ii)
    % NOTE: unrotated segment is constructed with y - r, not y + r, because
    % rotation must be on top half of segment (it's confusing)
    p = patch([x(ii); x(ii) + l(ii); x(ii) + l(ii); x(ii)],[y(ii); y(ii); y(ii) - r(ii); y(ii) - r(ii)], 'k', 'EdgeColor', clr, 'FaceColor', 'w', 'LineWidth', 2);
    
    % rotate given angle th(ii)
    rotate(p,[0 0 1],th(ii),[x(ii),y(ii),0]);
end
end