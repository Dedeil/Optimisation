function G = matGamma(siz,type)
%% matGamma     matrice Gamma de différentiation d'ordre 1 ou 2
%       G = matGamma(SIZE,TYPE)     crée la matrice de taille SIZE de   
%       dérivation de type TYPE
%       . SIZE = size(x) avec x le signal d'origine ou l'image d'origine non vectorisée
%       . TYPE = 'identity','gradient','laplacian'
%                                      
% 
%     Ex 1 (1D): 
%       t = [-10:.1:10]';
%       s = sin(t);
%       D = matGamma(size(s),'gradient');
% 
%       s_grad = D*s;
% 
%       figure(1); clf;
%       plot(t,s,'-',t,s_grad,'--');
%       legend('signal original','signal dérivé');
% 
% 
%     Ex 2 (2D):
%       im = im2double(imread('cameraman.tif'));
%       D  = matGamma(size(im),'gradient');
%       
%       im_grad  = D*im(:);
%       im_gradx = reshape(im_grad(1:length(im_grad)/2)    ,size(im));
%       im_grady = reshape(im_grad(length(im_grad)/2+1:end),size(im));
% 
%       figure(2); clf;
%       subplot(131); imshow(im,[]);       title('image originale');
%       subplot(132); imshow(im_gradx,[]); title('Gradient selon x de l''image');
%       subplot(133); imshow(im_grady,[]); title('Gradient selon y de l''image');
% 
% 
%     Ex 3 (2D):
%       im = im2double(imread('cameraman.tif'));
%       L  = matGamma(size(im),'laplacian');
% 
%       im_lap = reshape(L*im(:),size(im));
% 
%       figure(3); clf;
%       subplot(121); imshow(im,[]);     title('image originale');
%       subplot(122); imshow(im_lap,[]); title('Laplacien de l''image');
%%

dim = (min(siz) > 1) + 1;
if dim == 1, siz = sort(siz,'descend'); end

S = siz(1)*siz(2);

switch type
    case {'identity','Identity'}
        G = speye(S);
        
    case {'gradient','Gradient'}
        Dx = spdiags([-ones(S,1) ones(S,1)], [0,siz(1)], S, S)/2;
        Dx(end-siz(1)+1:end,:) = 0;
        
        Dy = spdiags([-ones(S,1) ones(S,1)], [0,1], S, S)/2;
        Dy(siz(1):siz(1):end,:) = 0;
        
        switch dim
            case 1, G = Dy;
            case 2, G = [Dx; Dy];
        end
        
    case {'laplacian','Laplacian'}  % V4
        a = repmat([0; ones(siz(1)-1,1)],[siz(2),1]);   a = [a(2:end); nan];
        b = repmat([ones(siz(1)-1,1); 0],[siz(2),1]);   b = [nan; b(1:end-1)];
        l = [zeros(siz(1),1); ones(S-siz(1),1)];        l = [l(siz(1)+1:end); nan(siz(1),1)];
        r = [ones(S-siz(1),1); zeros(siz(1),1)];        r = [nan(siz(1),1); r(1:end-siz(1))];
        
        Ly = spdiags([a b], [-1,1], S, S);
        Lx = spdiags([l r], [-siz(1),siz(1)], S, S);
        
        L = (dim==2)*Lx + Ly;
        L = L - spdiags(sum(L,2), 0, S, S);
        
        G = L;
        
    otherwise
        error('Type de différentiation inconnu.');
end