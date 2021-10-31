function [mainTiling, NFMAIN] = getMesoVoidPolygons(cx,cy,cij,L)
%% FUNCTION to construct void tilings from contact network

% mod cx and cy
cx = mod(cx,L);
cy = mod(cy,L);
NCELLS = length(cx);

%% Construct contact network across periodic images

% construct full list of cell centers
NP = 9*NCELLS;
P = zeros(NP,2);
for nn = 1:NCELLS
    cxtmp = cx(nn);
    cytmp = cy(nn);
    kbox = 1;
    for xx = -1:1
        for yy = -1:1
            P((kbox-1)*NCELLS + nn,1) = cxtmp + xx*L;
            P((kbox-1)*NCELLS + nn,2) = cytmp + yy*L;
            kbox = kbox + 1;
        end
    end
end

% loop over all possible image contacts
fprintf('\t ** Constructing adjacency matrix ...\n');
A = zeros(NP);
for nn = 1:NCELLS
    xn = cx(nn);
    yn = cy(nn);
    for mm = nn+1:NCELLS
        if cij(nn,mm) > 0
            xm = cx(mm);
            ym = cy(mm);

            % check if image contact
            dx = xm - xn;
            imx = round(dx/L);

            dy = ym - yn;
            imy = round(dy/L);
            
            % loop over image tiles, only add non-boundary contacts
            for xx = -1:1
                xi = 2 + xx;
                for yy = -1:1
                    yi = 2 + yy;
                    bi = 3*(xi-1) + yi;
                    
                    % add to periodic planar graph if not on exterior boundary
                    ytop = (imy == -1 && yy == 1);
                    xright = (imx == -1 && xx == 1);
                    ybottom = (imy == 1 && yy == -1);
                    xleft = (imx == 1 && xx == -1);
                    
                    if (~ytop && ~ybottom && ~xleft && ~xright)
                        bimm = 3*(xi - imx - 1) + yi - imy;
                        if (bimm < 0)
                            test = 1;
                        end
                        
                        nnInd = (bi-1)*NCELLS + nn;
                        mmInd = (bimm-1)*NCELLS + mm;
                        A(nnInd,mmInd) = 1;
                        A(mmInd,nnInd) = 1;
                    end
                end
            end
        end
    end
end

% create graph of whole periodic boundary
G = graph(A);

% get spatial graphs
fprintf('\t ** Getting spatial graph...\n');
obj = spatialgraph2D(G,P(:,1),P(:,2));
pgon = polyshape(obj);


% get polygons in main cell
fprintf('\t ** Reconstructing polygons in main cell...\n');
vpos = {pgon.Vertices}';
NF = size(vpos,1);
mainTiling = cell(NF,1);
NFMAIN = 0;
for ff = 1:NF
    vptmp = vpos{ff};
    cxtmp = mean(vptmp(:,1));
    cytmp = mean(vptmp(:,2));
    if cxtmp > 0 && cxtmp < L && cytmp > 0 && cytmp < L
        NFMAIN = NFMAIN + 1;
        mainTiling{NFMAIN} = vptmp;
    end
end
mainTiling(NFMAIN:end) = [];

fprintf('\t ** Finished!\n');

end