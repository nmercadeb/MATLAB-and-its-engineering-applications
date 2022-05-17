function [cos1,sin1,cos2,sin2,sin11,cos11,sin22,cos22] = angles(fill1,fill2,X,xf,n,u,v,vs)
%
% This function stablishes the direction of the daugther branches, having
% defined the mother's direction, the separation the bifurcation vectors.

% Flow in each branch
q1 = sum(X(:,4)==fill1);
q2 = sum(X(:,4)==fill2);
qT = q1 + q2;
q1 = q1/qT;
q2 = q2/qT;

% Values of alpha that minimazes the work
cos1=(1+q1^(4/n)-(1-q1)^(4/n))/(2*q1^(2/n));
cos2=(1+q2^(4/n)-(1-q2)^(4/n))/(2*q2^(2/n));

% Mass center of the irrigated volumes
cm1 = [0 0 0] + mean(X(X(:,4)==fill1,1:3));
cm2 = [0 0 0] + mean(X(X(:,4)==fill2,1:3));

% Projection of the mass center to the bifurcation plane
A=u(1); B=u(2); C=u(3); D=-(A*xf(1)+B*xf(2)+C*xf(3)); % Bifurcation plane

if A*cm1(1)+B*cm1(2)+C*cm1(3)+D==0 % If it is 0, cm is in the plane
    cmp1=cm1;
else
    cmp1=punt_pla(cm1,u,xf); % If not, compute the projection
end

if A*cm2(1)+B*cm2(2)+C*cm2(3)+D==0
    cmp2=cm2;
else
    cmp2=punt_pla(cm2,u,xf);
end

% Calculate the final angle within the bifurcation plane by calculating the mean
% Daugther 1
angle_alpha = acos(cos1); % Angle that minimizes the work
vect = cmp1-xf;
angle_beta = acos((vect*v')/(norm(vect)*norm(v)));  % Angle between mother and projecion of the mass center 
angle = min(mean([angle_alpha,angle_beta]),pi/2); % Angle of bifurcation within the bifrucation plane for daugther 1
cos1 = cos(angle);
sin1 = sin(angle);

% Same for daugther 2
angle_alpha = acos(cos2);
vect = cmp2-xf;
angle_beta = acos((vect*v')/(norm(vect)*norm(v)));
angle = min(0.5*angle_alpha+0.5*angle_beta,pi/2);
cos2 = cos(angle);
sin2 = sin(angle);


% New implementation: non coplanar bifurcations, branches can have 2
% degrees of freedom. Thus, the new branches can move in the separation
% plane folllwing the center of mass of the irrigated volume:

A=vs(1); B=vs(2); C=vs(3); D = -(A*xf(1)+B*xf(2)+C*xf(3)); % Separation plane

% Projection of the mass center into the separation plane
if A*cm1(1)+B*cm1(2)+C*cm1(3)+D==0 

    cmp1=cm1;
else
    cmp1=punt_pla(cm1,vs,xf);
end

if A*cm2(1)+B*cm2(2)+C*cm2(3)+D==0
    cmp2=cm2;
else
    cmp2=punt_pla(cm2,vs,xf);
end

% Finally the angle with respect to the mother vector:
vect = cmp1 - xf;
alpha = acos((vect*v')/(norm(vect)*norm(v)));
cos11 = cos(alpha); sin11 = sin(alpha);

vect = cmp2 - xf;
alpha = acos((vect*v')/(norm(vect)*norm(v)));
cos22 = cos(alpha); sin22 = sin(alpha);


end