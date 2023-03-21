function imf = opB(x)
%%

dim = (min(size(x)) > 1) + 1;


switch dim
    case 1
%         w = length(x);
%         U = -w/2+1/2:w/2-1/2;
%         dist  = sqrt(U.^2/w^2);      % anisotropic distance
%         n  = 5;
%         fc = 1/3;       % in [0,1] since D is normalized
%         Lo  = 1./(1 + (dist./fc).^(2*n));
%         
%         imf = ifft(ifftshift(Lo.*fftshift(fft(x))));
        error('Not implemented yet');
        
    case 2
        [h,w] = size(x);
        [U,V] = meshgrid(-w/2+1/2:w/2-1/2, -h/2+1/2:h/2-1/2);
        dist  = sqrt(U.^2/w^2 + V.^2/h^2);      % anisotropic distance
        n  = 5;
        fc = 1/3;       % in [0,1] since D is normalized
        Lo  = 1./(1 + (dist./fc).^(2*n));
        
        imf = ifft2(ifftshift(Lo.*fftshift(fft2(x))));
end




end