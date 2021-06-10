function curvaturesOfSheet = estimateCurvatureLS(referenceNodes,currentNodes,polynomToUse)
%function curvatureOfSheet(referenceNodes,currentNodes) estimates the
%curvature of the sheet strip by fitting the point cloud of nodal
%displacements in the z-direction as a function of the nodal position to a
%polynom:
%
%   w(x,y) = C_1*x^2 + C_2*x*y + C_3*y^2
%



% Renormalization
x1 = 1e-6*referenceNodes(:,2);
y1 = 1e-6*referenceNodes(:,3);
z2 = 1e-6*currentNodes(:,4);

x1 = x1 - min(x1);
y1 = y1 - min(y1);


if strcmp(polynomToUse,'2ndOrder-full')
    Xmat = [x1.^2 x1.*y1 y1.^2 x1 y1 ones(size(x1))];
    
elseif strcmp(polynomToUse,'2ndOrder-reduced')
    Xmat = [x1.^2 x1.*y1 y1.^2 ];
    
else
    disp('not implemented, see documentation.')
    
end

CParameters = lsqlin(Xmat,z2,[],[]);


if strcmp(polynomToUse,'2ndOrder-full')
    wFitted = CParameters(1)*x1.^2 + CParameters(2)*x1.*y1 + CParameters(3)*y1.^2 + CParameters(4)*x1 + CParameters(5)*y1 + CParameters(6)*ones(size(x1));
    wFittedFcn = @(a,b) CParameters(1)*a.^2 + CParameters(2)*a.*b + CParameters(3)*b.^2 + CParameters(4)*a + CParameters(5)*b + CParameters(6)*1;
elseif strcmp(polynomToUse,'2ndOrder-reduced')
    wFitted = CParameters(1)*x1.^2 + CParameters(2)*x1.*y1 + CParameters(3)*y1.^2 ;
    wFittedFcn = @(a,b) CParameters(1)*a.^2 + CParameters(2)*a.*b + CParameters(3)*b.^2 ;
else
    disp('not implemented, see documentation.')
    
end

[X,Y] = meshgrid(linspace(min(x1), max(x1),100),linspace(min(y1), max(y1),100));

kappaXX = -2 * CParameters(1);
kappaYY = -2 * CParameters(3);
kappaXY = - CParameters(2);

curvaturesOfSheet = [kappaXX kappaYY kappaXY];



% temp = yhatFcn(X,Y);
% figure;
% surf(X,Y,yhatFcn(X,Y))
% xlabel x; ylabel y; zlabel z
% hold on
% plot3(x1,y1,z2,'.b') 
% zlim(10e-4*[-1 1])

% folders = subdirImport('plots\','regex','R');
% print(horzcat('plots\R',num2str(length(folders))),'-dpng')
% close;








