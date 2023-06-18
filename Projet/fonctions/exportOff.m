function exportOff(filemane,coords,faces)
%% exportOff    exporte un maillage au format .off
%   exportOff(FILENAME,COORDS,FACES) exporte le maillage triangulé, dont 
%   les coordonnées 3D sont stockées dans la matrice Nx3 COORDS, et les
%   triangles dans la matrice Nx3 FACES, dans un fichier FILENAME (au format 
%   mon_fichier.off).
%%

f  = fopen(filemane,'w+');

% header
fprintf(f,'# Geometrical Mesh [NO_NAME] in [save_off.off]\n');
fprintf(f,'# Mesh Exported on %s\n',date);
fprintf(f,'#\n');
fprintf(f,'# N_vertices [%i]\n',size(coords,1));
fprintf(f,'# N_polygons [%i]\n',size(faces,1));
fprintf(f,'#\n');
fprintf(f,'OFF\n');
fprintf(f,'%i %i 0\n',size(coords,1),size(faces,1));

% print coords
dlmwrite(fn, coords, 'delimiter','\t','precision',5,'-append');
pause(.1);

% print faces
tri = [3*ones(size(faces,1),1), faces];
dlmwrite(fn, tri, 'delimiter','\t','precision',5,'-append');
pause(.1);

% footer
%fprintf(f,'EOF Mesh Geometry save_off.off [OK]\n');    % does not work together with dlmwrite 
dlmwrite(fn, 'EOF Mesh Geometry save_off.off [OK]','-append','delimiter','');

fclose(f);

end