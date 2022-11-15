function drawSSProjPoly(x, y, r, clr, Lx, Ly, pbcs)
% number of vertices
n = length(x);

% segments
ip1 = [2:n 1];
im1 = [n 1:n-1];
lx = x(ip1) - x;
ly = y(ip1) - y;
l = sqrt(lx.^2 + ly.^2);
nx = ly ./ l;
ny = -lx ./ l;

if pbcs == 1
    dbox = -1:1;
else
    dbox = 0;
end

% draw patch
for xx = dbox
    for yy = dbox
        xproj = x + r .* 0.5 .* (nx + nx(im1));
        yproj = y + r .* 0.5 .* (ny + ny(im1));
        patch('Faces',[1:n 1],'Vertices',[xproj + xx*Lx, yproj + yy*Ly],'FaceColor',clr,'EdgeColor','k');
%         plot(x, y, 'ko', 'markersize', 4, 'markerfacecolor', 'k');
    end
end

end