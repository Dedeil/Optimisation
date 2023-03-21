function G = matGamma3D(coords,faces,type)
%% matGamma_3D     matrice Gamma de differentiation d'ordre 1 ou 2 en 3D
%       G = matGamma_3D(COORDS,FACES,TYPE)  cree la matrice de derivation 
%       de type TYPE adaptee au maillage (XYZ,FACES)
%       . COORDS    matrice Nx3 dont la i-eme ligne correspond au vecteur
%                   de coordonnees (X,Y,Z) du i-eme point du maillage
%       . FACES     matrice Nx3 dont la i-eme ligne donne les  numeros des
%                   3 points composant un triangle du maillage
%       . TYPE = 'gradient','laplacian'
%                                      
%
%     Ex :
%       [C, F] = loadOff('teapot.off');
% 
%       L = matGamma3D(C,F,'laplacian');
% 
%       C_lap = C;
%       for it = 1:15
%           C_lap  = C_lap + 0.01*L*C_lap;
%       end
% 
%       figure(1); clf;
%       subplot(121); trisurf(F,C(:,1),C(:,2),C(:,3));              title('maillage original');
%       subplot(122); trisurf(F,C_lap(:,1),C_lap(:,2),C_lap(:,3));	title('Laplacien du maillage');
%%


choice_weights = 'no';  % 'exp'

fprintf('Compute smoothing operator Gamma... ');
N = max(faces(:));   % number of nodes

%% Connectivity
degree    = zeros(N,1);
neighbors = cell(N,1);
weights   = cell(N,1);

nb_edges = 0;  
edges    = [];
W        = [];

for v = 1 : N
    % find vertices' neighborhood
    [i,~] = find(faces==v);
    neighborhood = faces(i,:);
    neighborhood = unique(neighborhood(:));
    degree(v)    = length(neighborhood) - 1;   % -1 to uncount the point itself
    neighbors{v} = neighborhood(neighborhood ~= v)';
    
    % find edges
    nb_edges = nb_edges + nnz( v>0 & v<neighbors{v} );
    
    tmp   = neighbors{v}(v < neighbors{v})';
    edges = [edges; v*ones(length(tmp),1), tmp];
    
    % create weights matrix
    switch choice_weights
        case 'no'
            weights{v} = ones(1,degree(v));
            
        case 'exp'  
            vdist = coords(neighborhood,:)-repmat(coords(v,:),length(neighborhood),1);
            dist  = sqrt(sum(vdist.*vdist,2));
            clear vdist
            
            %%% CHECK THIS
            weights{v} = exp(-dist(neighborhood ~= v)');
            
        otherwise
            weights{v} = ones(1,degree(v));
    end
    W = [W; weights{v}(v < neighbors{v})'];
    
  
end
%max_nb_neighbors = max(degree);

 
%% Generate smoothing operator
M = size(edges,1);  %nb_edges;
switch type
    case 'gradient'     % computed as the graph incidence matrix
        I = sparse([1:M,1:M],[edges(:,1),edges(:,2)],[W,-W],M,N);
        
        G = I;
        
        
    case 'laplacian'
        adj = zeros(N);
%         for m = 1:M
%             adj(edges(m,1),edges(m,2)) = 1; %W(m);
%             adj(edges(m,2),edges(m,1)) = 1; %W(m);
%         end
        adj(sub2ind(size(adj),edges(:,1),edges(:,2))) = 1; %W(m);
        adj(sub2ind(size(adj),edges(:,2),edges(:,1))) = 1; %W(m);

        d = diag(degree);
        L = d - adj;
        
        G = sparse(L);
        
    otherwise
        error('Not implemented yet.');
end


prop.nb_vertices = N;
prop.nb_edges    = M;
prop.degree      = degree;
prop.neighbors   = neighbors;
prop.edges       = edges;

fprintf('done.\n');

end