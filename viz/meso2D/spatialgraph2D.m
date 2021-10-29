classdef spatialgraph2D

    properties (SetAccess=protected)
        G %node graph
        x %x coordinates
        y %y coordinates
        labels %node labels 
    end
    
    properties
        pruneType='basic';
    end
    
    properties (Hidden)
        tol
    end
    
    methods
        function  obj=set.pruneType(obj,str)
        
            obj.pruneType=validatestring(str,{'basic','legacy'});
            
        end
    end
    
    
    methods
        
        function obj = spatialgraph2D(G,x,y,labels,tol)
        % spatialgraph2D(G,x,y,labels)
        % spatialgraph2D(G,[],[])
        % spatialgraph2D(G) 
        %
        %IN:
        %
        %           G: a graph object
        %         x,y: vectors of x- and y-coordinate data for the node 
        %              locations. If omitted or if x=[],y=[] then default 
        %              coordinates will be assigned as when doing plot(G).             
        %      labels: Optional numeric or string vector of node labels.
        %              Defaults to consecutive integers.
            

            if nargin==1 ||isempty(x) && isempty(y)
              Hfig=figure('Visible','off');
              Hax=axes(Hfig);
              Hg=plot(Hax,G);
              x=Hg.XData;
              y=Hg.YData;
              delete(Hfig);
            end
            
            x=x(:); y=y(:);
            s=G.Edges{:,1}(:,1);
            t=G.Edges{:,1}(:,2);
            
            N=G.numnodes;
            lengthCheck=numel(x)==N && numel(y)==N;
            assert(lengthCheck,'Input vectors x and y must have lengths equal to the number of graph nodes.');
            
            
            w=sqrt((x(s)-x(t)).^2 + (y(s)-y(t)).^2);
            
            if nargin<4
                labels=1:G.numnodes;
            end
            
            if nargin<5 || strcmp(tol,'auto')
                D=pdist2([x,y],[x,y],'Euclidean','Smallest',2);
                obj.tol=min(D(2,:))/1000;
            end
            
            obj.G=graph(s,t,w);
            obj.x=x; obj.y=y;
            obj.labels=labels;
            
        end
        
        function [pgon,labelIndices]=polyshape(obj)
        %Find constituent polygons.
        %
        %  [pgon,nodeIDs]=polyshape(obj)
        %
        %
        %OUT:
        %
        %       pgon: A polyshape array containing the constituent
        %             polygons of the graph, when viewed as a mosaic.
        %
        %    nodeIDs: A cell array of corresponding node IDs of the
        %             nodes forming each polygon.
            

        
           if obj.pruneType=="legacy"
            obj=prune(obj);
           else
            obj=pruneBranches(obj);               
           end
           
            [G,x,y,lp,tol]=deal(obj.G,obj.x,obj.y,obj.labels,obj.tol);
            edgeMatrix=G.Edges{:,1};
            edgeList=1:size(edgeMatrix,1);
            
            
            %Initialize
            
            Vall=[x(:),y(:)];
            
%             [P,E]=polyshortest(obj,1);
%             pgon(1)=polyshape(Vall(P,:),'Simplify',true);
%             %Vcell={Vall(P,:)};

%             labelIndices{1}=lp(P);
            
            assigned=[];
            unassigned=setdiff(edgeList,assigned);

            k=0;
            
            %Iterate
            while ~isempty(unassigned)
                
                k=k+1;
                
                E0=unassigned(1);
                
                [P,E]=polyshortest(obj,E0);
                if ~isempty(P)
                    
                 pgon{k}=polyshape(Vall(P,:),'Simplify',true);
                 labelIndices{k}=lp(P);
                 assigned=union(assigned,E);

                else
                    assigned=union(assigned,E0);
                end
                   unassigned=setdiff(edgeList,assigned);       
            end
            pgon=[pgon{:}];
            
            
            %U=union(pgon);
            %U=union(regions(union(pgon)));
            %U=union(polybuffer(U,tol));
            
            U=union(polybuffer(pgon,tol));
            
            
            
            if ~U.NumHoles, return; end
            
            
            
            hgon=holes(U).';
            
            for i=1:numel(hgon)
                P= knnsearch(Vall,hgon(i).Vertices);
                pgon(end+1)=polyshape(Vall(P,:),'Simplify',true);
                labelIndices{end+1}=lp(P.'); %#ok<*AGROW>
            end
            
            
            
        end
        
        function varargout=plot(obj)
        %Display the graph. The nodes are positioned according to the
        %x and y data used to construct the object.
        %
        % Hg=plot(obj)
        %
        %OUT:
        %
        %  Hg: the PlotGraph handle
            
            h = plot(obj.G,'XData',obj.x,'YData',obj.y,'linewidth',2,'MarkerSize',7);
            nl = (obj.labels)+"";
            h.NodeLabel = '';
            xd = get(h, 'XData');
            yd = get(h, 'YData');
            text(xd, yd, nl, 'FontSize',17, 'FontWeight','bold', ...
                'HorizontalAlignment','left', 'VerticalAlignment','top');
            set(gca,'Fontsize',15,'FontWeight','Bold','LineWidth',2, 'box','on');
            
            if nargout, varargout={h}; end
            
        end
        
        function varargout=mosaic(obj)
        %Display the graph with constituent polygons shaded.
        %
        % [Hg,Hp]=mosaic(obj)
        %
        %OUT:
        %
        %  Hg: the PlotGraph handle
        %  Hp: Polygon graphics object
            
            Hg=plot(obj);
            
            pgon=polyshape(obj);
            hold on
            Hp=plot(pgon);
            hold off
            
            if nargout, varargout={Hg,Hp}; end
        end

    end
    
    methods (Hidden)
        
        function obj = prune(obj)
           
            switch obj.pruneType
               
                case 'basic'
                    
                    obj = pruneBasic(obj);
                    
                case 'legacy'
                
                    obj = pruneLegacy(obj);
                
            end
            
            
        end
        
        function obj = pruneBranches(obj)
        %Remove dangling branches
            
            Gp=obj.G;
            xp=obj.x;
            yp=obj.y;
            lp=obj.labels;

            
            tips=find(degree(Gp)<=1);
            
            while ~isempty(tips)
                
                Gp=rmnode(Gp,tips);
                xp(tips)=[];
                yp(tips)=[];
                lp(tips)=[];
                tips=find(degree(Gp)<=1);
                
            end

            
            obj=spatialgraph2D(Gp,xp,yp,lp);

        end
        
        function obj = pruneBasic(obj)
        %Remove dangling branches
            
        
            obj=pruneBranches(obj);
        
            Gp=obj.G;
            xp=obj.x;
            yp=obj.y;
            lp=obj.labels;
            

            N=Gp.numedges;
            onecuts=false(N,1);
            for i=1:N
                
                onecuts(i)=isempty( polyshortest(obj,i) );
                
            end
            
            Gp=Gp.rmedge(find(onecuts)); %#ok<FNDSB>
            
            
            tips=find(degree(Gp)<=1);

                Gp=rmnode(Gp,tips);
                xp(tips)=[];
                yp(tips)=[];
                lp(tips)=[];

            
            obj=spatialgraph2D(Gp,xp,yp,lp);

        end
        
        
        function obj = pruneLegacy(obj)
        %Remove dangling branches and polygons attached to them
 
        
            obj=pruneBranches(obj);
        
            Gp=obj.G;
            xp=obj.x;
            yp=obj.y;
            lp=obj.labels;
 
            bins=biconncomp(Gp);
            [~,cmaj]=max(histcounts(bins,1:max(bins)+1));
            ekeep=find(bins==cmaj);
            [sk,tk]=findedge(Gp,ekeep);
            [uk,~,jk]=unique([sk,tk]);
            xp=xp(uk); yp=yp(uk); lp=lp(uk);
            sk=jk(1:end/2); tk=jk(end/2+1:end);
            
            obj=spatialgraph2D(graph(sk,tk),xp,yp,lp);
            
            
        end
        
        function [P,E]=polyshortest(obj,arg1,arg2)
            
            G=obj.G;
            
            if nargin==3
                [nodeA,nodeB]=deal(arg1,arg2);
                %edgeNum=findedge(G,nodeA,nodeB);
            elseif nargin==2
                edgeNum=arg1;
                nodeA=G.Edges{edgeNum,1}(1);
                nodeB=G.Edges{edgeNum,1}(2);
            else
                error 'Too few arguments'
            end
            
            Gtmp=rmedge(G,nodeA,nodeB);
            P=shortestpath(Gtmp,nodeA,nodeB);
            
            if isempty(P), E=[]; return ;end
            
            E=G.findedge(P,[P(2:end),P(1)]);
        end
        
    end
    
    methods %"Inherited" from graph objects
        
        
    end
    
end
