function edgeidx = isedge(bij,bw)
    NVTOT = length(bij);
    [n,~,N] = size(bw);
    edgeidx = zeros(NVTOT,1);
    gi = 1;
    for cc = 1:N
        for vv = 1:n
            % see if gi and gip1 are connected to different things
            
            
    
    
    for nn = 1:N
        for ii = 1:n
            if bw(ii,1,nn) ~= 0 && edgeidx(gi) == 1
                edgeidx(gi) = 0;
            end
            gi = gi + 1;
        end
    end
end