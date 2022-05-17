% Report figures:

% 1. Separation plane within the right lung between volumes 1:3 and 4:11:
load xarxa.mat

Xn = Xarxa;
Xn(ismember(SE,1:3),4)  = 3;
Xn(ismember(SE,4:11),4) = 4;

g1 = Xn(:,4) == 3; % Segment 1
g2 = Xn(:,4) == 4; % Segment 2

% Separation plane computation
X = [Xn(g1,1:3);Xn(g2,1:3)];
seg = cell(length(X),1); 
seg(1:sum(g1)) = {'s1'}; seg((1:sum(g2))+ sum(g1)) = {'s2'};
MdlLinear = fitcdiscr(X,seg);
% Plane equation is: Ax + By + Cz + D = 0
K = MdlLinear.Coeffs(1,2).Const;  % D
L = MdlLinear.Coeffs(1,2).Linear; % Separation plane vector (non unitary): [A, B, C]
vs = L/norm(L);

% Representation
scatter3(Xn(g1,1),Xn(g1,2),Xn(g1,3)); hold on % Segment 1
scatter3(Xn(g2,1),Xn(g2,2),Xn(g2,3));         % Segment 2 
xlabel('x');ylabel('y');zlabel('z');
f = @(x1,x2) (-1/L(3))*(K + L(1)*x1 + L(2)*x2); % Plane equation as z = -1/(C*(Ax+By))
h = fsurf(f,'FaceColor','k');               % Plane surface in black





