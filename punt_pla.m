function projeccio = punt_pla(punt, normal, puntdelpla)
%
% This function projects the point 'punt' within the plane defined by the
% orthogonal vector 'normal' and the point 'puntpla'. The output are the
% coordinates of the point within the plane.

A=normal(1); B=normal(2); C=normal(3);
D=-(A*puntdelpla(1)+B*puntdelpla(2)+C*puntdelpla(3));

t=-(A*punt(1)+B*punt(2)+C*punt(3)+D)/(A*A+B*B+C*C);

projeccio(1)=punt(1)+A*t;
projeccio(2)=punt(2)+B*t;
projeccio(3)=punt(3)+C*t;

end