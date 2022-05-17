clear all, close all

% First we load data from the lungs:
%   1. segments: is a 512x512x251 matrix filled with integers from 0 to 20
%      defining each lobe of the lung (1-->20) and the outside of the lung
%      (0)
%   2. Xarxa: is a matrix obtained from the previous one and addition data.
%      The first three columns of the matrix contain cartesina coordinates
%      for each point of 2, this has been computed with addition data that
%      stablished the ditance in milimeteres between points in 1, in the x,
%      y and z directions. The fourth column of this matrix states numbers
%      1 or 2 if it correspond to the left or right lung, and the fifth
%      column has numbers from 1 to 20 for each lobe.    

load segments.mat % 1
load xarxa.mat                                        % 2

% Separate the matrix Xarxa into a the column stablishing the segments and 
% the another matrix containing the cooridnates for each point and a fourth
% column to identify the points.
SE = Xarxa(:,5);      % Segments numbered (1 to 20)
Xarxa = Xarxa(:,1:4); % Forth colum 1 or 2 for right and left lung respectively

% General parameters
M = 3; % Proportional factor for diameter-length of the branload ("Segments_Macacs_SenseFiltres.mat",'segments');ch
n = 3; % Optimum factor for the function Angles

% TRACHEA
dTraq = 7.8; % Tracha diameter

% Matrix contaning all branches: each row is a bronchial branch
% Meaning of each column:
%  1: Howe many branches there are between this one and the trachea
%  2: In which row is its mother branch
%  3: Diameter
%  4: Longitude
%  5: xi, initial x coordinate
%  6: yi, initial y coordinate
%  7: zi, initial z coordinate
%  8: xf, final x coordinate
%  9: yf, final y coordinate
% 10: zf, final z coordinate
% 11: birth, number of daughter branches emerging from it
% 12: vx, x coordinate of its mother's vector
% 13: vy, y coordinate of its mother's vector
% 14: vz, z coordinate of its mother's vector
% 15: angle, angle of bifrucation of its mother's branch
Tubs = [0 0 dTraq 3*dTraq 0 0 3*dTraq 0 0 0 2 -1 0 0 0];

% FIRST BIFRUCATION -------------------------------------------------------
N1 = sum(Xarxa(:,4)==1); % Volume right lung
N0 = sum(Xarxa(:,4)==1) + sum(Xarxa(:,4)==2); % Total volume

d1 = dTraq*(N1/N0)^(1/3); % Diameter for the branch irrigating the right lung
d2 = dTraq*(sum(Xarxa(:,4)==2)/size(Xarxa,1))^(1/3); % Diameter for the left one

l1 = M*d1; % Longitud of the right branch
l2 = M*d2; % Left branch

% fill1 and fill2 define the value in hte fourth colume of Xarxa that
% reference the points of the irrigated volume for each branch
fill1 = 1; 
fill2 = 2;
X     = Xarxa;
xf    = [0 0 0];  % End point of the mother branch
u     = [0 -1 0]; % Unitary orthogonal vector to the bifurcation plane
v     = [0 0 -1]; % Unitary vector in the directon of the mother's branch
vs    = [-1 0 0]; % Unitary orthogonal vector to the separation plane
[cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22] = angles(fill1,fill2,X,xf,n,u,v,vs);
% cos and sin with one number reference the direction in the bifrucation
% plane, following the axes defined by the mother branch and the separation
% plane. The ones with two numbers define the direction in the separation
% plane, with respect to the axes defined by the mother branch and the
% bifurcation plane.

dir1 = v*(cos1+cos11)+vs*sin1-u*sin11; dir1 = dir1/norm(dir1); % direcction daugther 1
dir2 = v*(cos2+cos22)-vs*sin2-u*sin22; dir2 = dir2/norm(dir2); % direcction daugther 2

% The new branches are:    
Tubs = [Tubs;
    1 1 d1 l1 0 0 0 (v*cos1+vs*sin1)*l1 2 v acosd(cos1);
    1 1 d2 l2 0 0 0 (v*cos2-vs*sin2)*l2 2 v acosd(cos2)];

% RIGHT LUNG --------------------------------------------------------------

% Each colum of the cell array contains the referenced lobuls of each
% partition in each iteration. For example, in the first iteration will
% compute the branches that irrigate volume 1:3 and 4:11. The following
% iteration starts from the branch that irrigates the volume defined in the
% first row, in the case mentioned, it is 1:3 and is devided into 1:2 and 3
sec = { 1:3, 1:2, 1,  4:5, 4, 6:10, 7:10,  7:8, 7, 9;
       4:11,   3, 2, 6:11, 5,   11,    6, 9:10, 8, 10}; 

for k = 1:size(sec,2)
    Xn = Xarxa;
    Xn(ismember(SE,sec{1,k}),4) = 3; % Volume irrigated by daugther branch 1 is referenced with 3
    Xn(ismember(SE,sec{2,k}),4) = 4; % Volume irrigated by daugther branch 2 is referenced with 4

    g1 = Xn(:,4) == 3; % Logical indexes for volume 1
    g2 = Xn(:,4) == 4; % Logical indexes for volume 2

    % Separation plane is find by Linear Discriminant Analisis method
    X = [Xn(g1,1:3); Xn(g2,1:3)];
    seg = cell(length(X),1);
    seg(1:sum(g1)) = {'s1'}; % Points corresponding to class 1 (volume 1)
    seg((1:sum(g2))+ sum(g1)) = {'s2'}; % Points corresponding to class 2 (volume 2)
    MdlLinear = fitcdiscr(X,seg); % Adjunt linear discriminant model
    L = MdlLinear.Coeffs(1,2).Linear; % Separation plane perpendicular vector
    vs = L'./norm(L); % Unitary normal separation plane vector   
    
    % For the angles function:
    fill1 = 3;
    fill2 = 4;
    X     = Xn;
    if ismember(k, [4 6 10]) % In this columns of sec there is a change in the bronchial path
        switch k
            case 4
                in = 5; % Row of the branch leading to 4:5
                d1 = Tubs(in,3); % Mother's diameter
            case 6
                in = 11; % Row of the branch leading to 6:11
                d1 = Tubs(in,3);
            case 10
                in = 19; % Row of the branch leading to 9:10
                d1 = Tubs(in,3);
        end
        xf = Tubs(in,8:10);
        v  = Tubs(in,8:10) - Tubs(in,5:7);    v = v/norm(v); % Unitary mother branch vector 
    else
        xf = Tubs(end-1,8:10);
        v  = Tubs(end-1,8:10) - Tubs(end-1,5:7);    v = v/norm(v); % Unitary mother branch vector
    end
    u     = cross(v,vs);                           u = u./norm(u); % Unitary normal bifurcation plane vector
    % Since separation vector does not necessarly be orhotgonal to mother 
    % and bifrucation vectors, the cross product of these two is performed
    % to obtain the perpedindicular vector:
    us = cross(u,v); % Unitary normal separation vector for the axis
    [cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22] = angles(fill1,fill2,X,xf,n,u,v,us);
    
    dir1 = v*(cos1+cos11)+us*sin1-u*sin11; dir1 = dir1/norm(dir1); % Direction daugther 1
    dir2 = v*(cos2+cos22)-us*sin2-u*sin22; dir2 = dir2/norm(dir2); % Direction daugther 2
    
    if ismember(k,1:3)
        g = 1+k;
        m = 2*k;
    elseif k == 4
        g = 3;
        m = 5;
    elseif ismember(k,5:6)
        g = 4;
        m = k + 5;
    elseif k == 7
        g = 5;
        m = 14;
    elseif ismember(k,8:10)
        m = k + 9;
        if k == 8
            g = 6;
        else
            g = 7;
        end
    end

     % Volums for each branch = sum number of points
    N1 = sum(Xn(:,4)==3);
    N2 = sum(Xn(:,4)==4);
    N0 = N1 + N2;
    % Diametres
    d3 = d1*(N1/N0)^(1/3);
    d4 = d1*(N2/N0)^(1/3);    
    % Longituds
    l1 = M*d3;
    l2 = M*d4;

    d1 = d3; % For the next iteration

    Tubs = [Tubs;
        g m d3 l1 xf xf+dir1*l1 2 v acosd(cos1);
        g m d4 l2 xf xf+dir2*l2 2 v acosd(cos2)];
end

% Change natality of the branches irrigating the lobes:
Tubs([7:9, 12, 13, 15, 17, 20:23],11) = 0; 

% REPRESENTATION ----------------------------------------------------------
pl = 1:11;  % Right lobes
ll = length(pl);
colormap jet;
cmap = colormap;  % Colors for the plot

figure(1);
title('Right lung bronchial tree','FontSize',15,'Interpreter','latex')
hold on, axis equal
i = 1;
for k = 1:size(Tubs,1)  
    % All branches faded and black, except those irrigating a lobe. These
    % will be opaque and have the color of the lobe which they irrigate.
    if Tubs(k,11) == 0
        branca(Tubs(k,5:7),Tubs(k,8:10), Tubs(k,3)/2, 100, cmap(i*floor(length(cmap)/ll),:), 1)
        i = i + 1;
    else
        branca(Tubs(k,5:7),Tubs(k,8:10), Tubs(k,3)/2, 100, 'k', .25)
    end
end

for k = 1:ll
    plcolor=cmap(k*floor(length(cmap)/ll),:); % Divide the colormap into as many colors as lobes representing
    R = double(ismember(segments,pl(k))); % Get the desired volume of the matrix segments
    fR = isosurface(R, 0.99);  % Create the patch object, isovalue = 0.99 since we want the contour
    fR.faces = fliplr(fR.faces); % Ensure normals point OUT
    % data extract from TC
    coordenades_carina = [253 279 171]; 
    dX = 0.310547; % mida x pixel en mm
    dY = 0.310547; % mida y pixel en mm
    dZ = 0.625; % mida z pixel en mm
    fR = transformar_coordenades(fR,coordenades_carina,[dX dY dZ]); % Escalate to real values
    patch(fR,'FaceColor',plcolor,'FaceAlpha',0.15,'EdgeColor','none') % Represent patch object
end


% Function to represent the volumes ---------------------------------------
function fv = transformar_coordenades(fv,xc,dx)
% This function escaltes a patch object fv 
fv.vertices = (fv.vertices-xc).*dx;
end