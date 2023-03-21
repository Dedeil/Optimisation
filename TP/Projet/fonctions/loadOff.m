function [coords, faces] = loadOff(filename)
% loadOff      charge un fichier .off
%   [COORDS,FACES] = loadOff(FILENAME) charge le fichier FILENAME
%   (au format mon_fichier.off) et en extrait les informations 3D :
%
%   Output:
%   . COORDS    matrice Nx3 dont la i-eme ligne correspond au vecteur de 
%               coordonnees (X,Y,Z) du i-eme point du maillage
%   . FACES     matrice Nx3 dont la i-eme ligne donne les  numeros des 3
%               points composant un triangle du maillage

    fprintf('Load data... ');
    fid1 = fopen(filename,'r'); 

    l = fgetl(fid1);
    while strcmp(l(1),'#')   
        l = fgetl(fid1); 
    end

    % get #vertices, #faces and #voxels
    l = fgetl(fid1); 
    n = str2num(l);
    V = n(1);
    F = n(2);

    % get coordinates
    coords = zeros(V,3);
    for itv = 1:V
        l       = fgetl(fid1);
        coords(itv,:) = str2num(l);
    end

    % get faces
    faces = zeros(F,4);
    for itf = 1:F
        l         = fgetl(fid1);
        faces(itf,:) = str2num(l);
    end
    faces = faces(:,2:end)+1;   % faces needs to be >0 in triangulation.m
                                % Rq: 1st col gives the nb of vertices per face

    fprintf('done.\n');                            
end
