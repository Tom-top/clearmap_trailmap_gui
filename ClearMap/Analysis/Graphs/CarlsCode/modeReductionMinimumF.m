% mode reduction : add vectors one by one starting from the smallest and
% check for linear independence
% needs: findBasisVectors_fastF, bfs, isMatrixSingular2F


function modeReductionMinimumF(whichLattice,latticeSize)
clf
user_name =  java.lang.System.getProperty('user.name');
whichUser = char(user_name);

%whichLattice = 'RandomPoints3D'; %'Torus2holes';%'Torus1hole';%'RandomFloors3D';
%latticeSize = '1034';


mainfolder = ['/Users/' whichUser '/dropbox/currently_working_on/ModesProject'];

path(path,[mainfolder '/matlabmodesfunctions'])
path(path,[mainfolder '/2.findbasisvectors'])
%addpath(genpath([mainfolder '/matlabmodesfunctions']))
path(path,[mainfolder '/1.maketoygraphs'])
path(path,[mainfolder '/1.maketoygraphs/dijkstra'])
path(path,[mainfolder '/1.maketoygraphs/makeToyGraphsF/'])
path(path,[mainfolder '/2.findbasisvectors/findBasisVectorsF/'])
path(path,[mainfolder '/matlabmodesfunctions/matlab_bgl'])

saveinfolder = [mainfolder '/1.makeToygraphs/test' whichLattice '/'];
cd(mainfolder)

fileRead = ['test' whichLattice '_'  latticeSize  ]; %'T1'


% load data
load([saveinfolder  fileRead '.mat'],'LEAF');


% displayNetwork3DF(LEAF,1,@(x) (x/max(x)).^0.3)
% axis equal
% hold on

NVertex = length(LEAF.Vertex);
NBonds = length(LEAF.Bond);

% find vectors from all spaning trees rooted at each and every node 
 cellBasis = cell(1,NVertex);
 cellBasisO = cell(1,NVertex);
 %cellBasissum= zeros(1,NVertex);
tic

for jnode = 1:NVertex
    AdjacencyMat = sparse(~~LEAF.Conductivity);
    [d dt pred] = bfs(AdjacencyMat,jnode);
    pred(jnode)=[];
    Knodes = [(1:jnode-1) (jnode+1:NVertex)];
    spanningTree = sparse([Knodes pred], [pred Knodes], ones(1,2*NVertex-2));
    basisG = AdjacencyMat - spanningTree;
    
    [Krows Kcol] = find(triu(basisG > 0));
    
    
    
    Ebasis = findBasisVectors_fastF(spanningTree, Kcol, Krows, LEAF);
    cellBasisO{jnode} = sparse(Ebasis);
    %Ebasis = rotateBasis_Horton1F(Ebasis,LEAF, 1, 150 ,'n');
    disp([int2str(jnode) ' ' int2str(sum(Ebasis(:)))])
    cellBasis{jnode} = sparse(Ebasis);
    %cellBasissum(jnode) = sum(Ebasis(:));
end
toc

%[minfromTree imin] = min(cellBasissum);

EallvectorsO = cell2mat(cellBasisO);
EallvectorsO = unique(EallvectorsO','rows')'; % all the vectors
% rows are the basis vectors (matlab: (rows, cols))


EbasisMin =[];
VectorsF = [];
lengthEall = full(sum(EallvectorsO));


% check for linear independence and add vectors one by one 
jl=3; % smallest cycle is 3 edges long
tic
while jl<=max(lengthEall) && size(EbasisMin,2)<NBonds-NVertex+1
    
    Ksize = lengthEall==jl;
    Vectors = EallvectorsO(:,Ksize);
    for jadd=1:sum(Ksize)
        if mod(jadd,20)==0, 
            [jl jadd sum(Ksize)]
        end
        VectorsTry = [VectorsF Vectors(:,jadd)];
        
        [isSingular VectorsTryO] = isMatrixSingular2F(full(VectorsTry));
        if isSingular==0
            VectorsF = VectorsTryO;
            EbasisMin = [EbasisMin Vectors(:,jadd)];
        end
    end
    disp(['iteration: ' int2str(jl) ' length basis: ' int2str(size(EbasisMin,2))])
    jl=jl+1;
end        
toc

%     VectorsF = findLinearlyIndepF([EbasisMin Vectors],'first');
%     sum(full(VectorsF))
%     VectorsF = findLinearlyIndepF(VectorsF,'first');
%     sum(full(VectorsF))
%     VectorsF = fliplr(findLinearlyIndepF(fliplr(VectorsF),'last'));
%     sum(full(VectorsF))
%     EbasisMin = VectorsF;
%     jin = jin+size(Vectors,2);
    
%end
% %%
% isSingularAll = isMatrixSingularF(EbasisMin(:,1:end)');
% isSingularWithoutLast = isMatrixSingularF(EbasisMin(:,1:end-1)');
% if isSingularAll
%     EbasisMin(:,end) = [];
% end

fileWrite = [fileRead 'MinVectors'];

save([saveinfolder '/' fileWrite '.mat'],'LEAF','cellBasisO','EbasisMin')






function [d dt pred] = bfs(A,u,target)
% BFS Compute breadth first search distances, times, and tree for a graph
%
% [d dt pred] = bfs(A,u) returns the distance (d) and the discover time
% (dt) for each vertex in the graph in a breadth first search 
% starting from vertex u.
%   d = dt(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% [...] = bfs(A,u,v) stops the bfs when it hits the vertex v
%
% Example:
%   load_gaimc_graph('bfs_example.mat') % use the dfs example from Boost
%   d = bfs(A,1)
%
% See also DFS

% David F. Gleich
% Copyright, Stanford University, 2008-20098

% History
% 2008-04-13: Initial coding

if ~exist('target','var') || isempty(full), target=0; end

if isstruct(A), rp=A.rp; ci=A.ci; 
else [rp ci]=sparse_to_csr(A); 
end

n=length(rp)-1; 
d=-1*ones(n,1); dt=-1*ones(n,1); pred=zeros(1,n);
sq=zeros(n,1); sqt=0; sqh=0; % search queue and search queue tail/head

% start bfs at u
sqt=sqt+1; sq(sqt)=u; 
t=0;  
d(u)=0; dt(u)=t; t=t+1; pred(u)=u;
while sqt-sqh>0
    sqh=sqh+1; v=sq(sqh); % pop v off the head of the queue
    for ri=rp(v):rp(v+1)-1
        w=ci(ri);
        if d(w)<0
            sqt=sqt+1; sq(sqt)=w; 
            d(w)=d(v)+1; dt(w)=t; t=t+1; pred(w)=v; 
            if w==target, return; end
        end
    end
end



function [rp ci ai ncol]=sparse_to_csr(A,varargin)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.
%
% [rp ci ai] = sparse_to_csr(i,j,v,n) returns a csr representation of the
% index sets i,j,v with n rows.
%
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A)
%
% See also CSR_TO_SPARSE, SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-07: Initial version
% 2008-04-24: Added triple array input
% 2009-05-01: Added ncol output
% 2009-05-15: Fixed triplet input

error(nargchk(1, 5, nargin, 'struct'))
retc = nargout>1; reta = nargout>2;

if nargin>1
    if nargin>4, ncol = varargin{4}; end
    nzi = A; nzj = varargin{1};
    if reta && length(varargin) > 2, nzv = varargin{2}; end    
    if nargin<4, n=max(nzi); else n=varargin{3}; end
    nz = length(A);
    if length(nzi) ~= length(nzj), error('gaimc:invalidInput',...
            'length of nzi (%i) not equal to length of nzj (%i)', nz, ...
            length(nzj)); 
    end
    if reta && length(varargin) < 3, error('gaimc:invalidInput',...
            'no value array passed for triplet input, see usage'); 
    end
    if ~isscalar(n), error('gaimc:invalidInput',...
            ['the 4th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
    if nargin < 5, ncol = max(nzj); 
    elseif ~isscalar(ncol), error('gaimc:invalidInput',...
            ['the 5th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
else
    n = size(A,1); nz = nnz(A); ncol = size(A,2);
    retc = nargout>1; reta = nargout>2;
    if reta,     [nzi nzj nzv] = find(A); 
    else         [nzi nzj] = find(A);
    end
end
if retc, ci = zeros(nz,1); end
if reta, ai = zeros(nz,1); end
rp = zeros(n+1,1);
for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
if ~retc && ~reta, rp=rp+1; return; end
for i=1:nz
    if reta, ai(rp(nzi(i))+1)=nzv(i); end
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;