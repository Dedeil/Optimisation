function H = matH(siz,type,varargin)
%% matH     matrice H de flou
%     H = matH(SIZE,TYPE,N)     crée la matrice de taille SIZE correspondant
%     au masque d'un flou de type TYPE et de taille N impaire
%     . SIZE = size(x) avec x le signal ou l'image d'origine non vectorisée
%     . TYPE = 'gaussian','uniform'
% 
%     H = matH(SIZE,'arbitrary',F)      crée la matrice de taille SIZE   
%     correspondant au masque F
%     . SIZE = size(x) avec x le signal ou l'image d'origine non vectorisée
%     . F est un vecteur (si x est un signal) ou une matrice (si x est une 
%       image), de taille impaire
% 
%     Ex 1 (1D): 
%       t = [-10:.1:10]';
%       s = sin(t);
%       H = matH(size(s),'uniform',25);
% 
%       s_blur = H*s;
% 
%       figure(1); clf;
%       plot(t,s,'-',t,s_blur,'--');
%       legend('signal original','signal flouté');
% 
% 
%     Ex 2 (1D): 
%       t = [-10:.1:10]';
%       s = sin(t);
%       H = matH(size(s),'arbitrary',[1:25]');
% 
%       figure(2); clf;
%       s_blur = H*s;
% 
%       plot(t,s,'-',t,s_blur,'--');
%       legend('signal original','signal flouté');
% 
% 
%     Ex 3 (2D):
%       im = im2double(imread('cameraman.tif'));
%       H  = matH(size(im),'gaussian',7);
% 
%       im_blur = reshape(H*im(:),size(im));
% 
%       figure(3); clf;
%       subplot(121); imshow(im,[]);      title('image originale');
%       subplot(122); imshow(im_blur,[]); title('image floutée');
%%

dim = (min(siz) > 1) + 1;

% check parameters
if dim == 1, siz = sort(siz,'descend'); end

if nargin ~= 3
    error('Nombre d''entrées incorrect, la fonction attend 3 arguments.'); 
else 
    if max(strcmp(type,{'Arbitrary','arbitrary'})) == 1
        N = size(varargin{1});
    else
        N = varargin{1};
        if numel(N) > 1
            error('La taille du masque doit être un entier avec les options ''gaussian'' et ''uniform''.');
        end
        
        if dim == 1, N = [N,1]; end
        if dim == 2, N = [N,N]; end
    end
end

if min(mod(N,2)) == 0, error('Les dimensions du masque doivent être impaires.'); end

if sum(N > siz) > 0, error('Le masque est plus grand que le signal/l''image.'); end


% 

    
S    = siz(1)*siz(2);
c    = (N+1)/2;
bnd  = (N-1)/2;

supp.i = -bnd(1):bnd(1);
supp.j = -bnd(2):bnd(2);

% mask
switch type
    case {'Gaussian','gaussian'}
        std = N(1)/6;
        
        if dim == 1
            x   = supp.i';
            arg = -(x.*x)/(2*std*std);
            
            F = exp(arg);
        else 
            [x,y] = meshgrid(supp.j,supp.i);
            arg   = -(x.*x + y.*y)/(2*std*std);
            
            F = exp(arg);
        end
        
        
    case {'Uniform','uniform'}
        F = ones(N);
        
    case {'Arbitrary','arbitrary'}
        F = varargin{1};
      
    otherwise
        error('Type de flou inconnu.');
end
        

% (sparse) matrix from mask
dnum = repmat(supp.i',[1,N(2)]) + repmat(supp.j*siz(1),[N(1),1]);

d = [];

for j = supp.j
    for i = supp.i
        
        neigh = [zeros(-j*siz(1),1) ; ...
                 repmat([zeros(-i,1); ones(siz(1)-abs(i),1); zeros(i,1)],[siz(2)-abs(j) 1]) ; ...
                 zeros(j*siz(1),1)];
        neigh = circshift(neigh,i+j*siz(1));
        
        d = [d, F(c(1)+i,c(2)+j)*neigh];
    end
end

H = spdiags(d,dnum(:),S,S);

% normalize filter
H = bsxfun(@rdivide, H, sum(H,2));

end