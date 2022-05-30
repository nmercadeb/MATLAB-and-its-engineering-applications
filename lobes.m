% This code generates the file lobes.mat that contains the patch objects to
% represent the 20 lung lobes.

load ("Segments_Macacs_SenseFiltres.mat",'segments');
ceg = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', ...
    'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20'};

for k = 1:20
    R3 = double(ismember(segments,k));
    fR3 = isosurface(R3, 0.99);  % Create the patch object
    fR3.faces = fliplr(fR3.faces);% Ensure normals point OUT
    coordenades_carina = [253 279 171]; %carina = [252 280 163];
    dX = 0.310547; % mida x pixel en mm
    dY = 0.310547; % mida y pixel en mm
    dZ = 0.625; % mida z pixel en mm
    eval([ceg{k} ' = transformar_coordenades(fR3,coordenades_carina,[dX dY dZ]);']);
end

save('lobs.mat','S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', ...
    'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20','-mat')


function fv = transformar_coordenades(fv,xc,dx)
fv.vertices = (fv.vertices-xc).*dx;
end
