function Z = opHt(x,type,varargin)
%% Hadj     operateur adjoint du flou
%     Z = Hadj(X,TYPE,N)  calcule la convolution X par l'adjoint d'un flou 
%     de type TYPE et de taille  N impaire
%     . TYPE = 'gaussian','uniform'
%
%     Autrement dit, Hadj(X,TYPE,N) calcule le produit matriciel :   
%           H'*x 
% 
% 
%     Ex 1 (2D):
%       im = im2double(imread('cameraman.tif'));
%       
%       ker = 'gaussian';  
%       z = opHadj(im,ker);
% 
%       figure(1); clf;
%       subplot(121); imshow(im,[]);      title('image originale');
%       subplot(122); imshow(z,[]);       title('image flout�e');
%%

siz = size(x);

dim = (min(size(x)) > 1) + 1;

% check parameters
if dim == 1, siz = sort(siz,'descend'); end

if nargin ~= 3
    error('Nombre d''entr�es incorrect, la fonction attend 3 arguments.'); 
else 
    if max(strcmp(type,{'Arbitrary','arbitrary'})) == 1
        N = size(varargin{1});
    else
        N = varargin{1};
        if numel(N) > 1
            error('La taille du masque doit �tre un entier avec les options ''gaussian'' et ''uniform''.');
        end
        
        if dim == 1, N = [N,1]; end
        if dim == 2, N = [N,N]; end
    end
end

if min(mod(N,2)) == 0, error('Les dimensions du masque doivent �tre impaires.'); end

if sum(N > siz) > 0, error('Le masque est plus grand que le signal/l''image.'); end


switch type
    case 'gaussian', ker = fspecial('gaussian',N,N(1)/6);  
    case 'uniform',  ker = ones(N);
    case 'none',     ker = 1;
end

ker = ker./sum(ker(:));

switch dim
    case 1
        if isrow(x), ker = ker'; end
        Z = real(ifft(conj(psf2otf(ker,size(x))).*fft(x)));
        
    case 2
        %Z = real(ifft2(conj(psf2otf(ker,size(x))).*fft2(x)));
        Z = imfilter(x,ker','replicate','same');
end

end
