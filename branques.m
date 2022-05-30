function branques(inici, final, radi, punts, color, transparencia)
%
% This function represents a bronchial branch as a cilinder given the input
% arguments:
%  - Initial point
%  - Final point
%  - Radi of the cilinder
%  - Number of points to compute the cilinder (if it is samll it has a polygonal shape)
%  - Color
%  - Alpha face (escalar between 0 and 1)

[X, Y, Z] = cylinder(radi,punts); % Default cylinder given by matlab of specified radi and number of points

for k = 1:size(inici,1)

    normal = final(k,:)-inici(k,:);    % Direction of the branch = cylinder axis
    longitud = norm(normal); % Longitude
    Zr = longitud*Z;                   % Adjust cylinder length to the desired one

    % The cylinder is oriented algong the Z direction [0 0 1], we want to
    % rotate this orientation algong the direction stablished by 'normal':
    zDirection = normal/norm(normal); % Branch axis

    % Create orthonormal basis:
    if abs(dot(zDirection, [0 1 0])) < 1
        xDirection = cross([0 1 0], zDirection);
        xDirection = xDirection/norm(xDirection);
    else
        xDirection = [1, 0, 0];
    end

    yDirection = cross(zDirection, xDirection);

    % Rotation matrix
    a = [xDirection',yDirection',zDirection'];

    % Aply rotation + translation to the coordinates given be @cylinder matlab
    % function:
    Xr = a(1,1)*X+a(1,2)*Y+a(1,3)*Zr+inici(k,1);
    Yr = a(2,1)*X+a(2,2)*Y+a(2,3)*Zr+inici(k,2);
    Zr = a(3,1)*X+a(3,2)*Y+a(3,3)*Zr+inici(k,3);


    h1 = surf(Xr,Yr,Zr); % Cylinder surface
    hold on
    set(h1,'FaceColor',color, 'EdgeColor','None');
    alpha(h1,transparencia);

    color_fosc=get(h1,'FaceColor')/3; % Perimeter of the bases less faded

    h2=plot3(Xr(1,:),Yr(1,:),Zr(1,:),'linewidth',1); % base 1 cilindre
    set(h2,'Color',color_fosc);
    h3=plot3(Xr(2,:),Yr(2,:),Zr(2,:),'linewidth',1); % base 2 cilindre
    set(h3,'Color',color_fosc);
end
end