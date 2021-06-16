function plotDensityEstimate(nodalData,custom,xDim,yDim)
tic
figure('name','densityPlot');
hh= scatterhist(nodalData(1:3:end,2),nodalData(1:3:end,3),'Kernel','on','marker','.','color',custom.Color(3,:));
xlabel('x [m]'); ylabel('y [m]');
hold on
plot([xDim(1) xDim(2) xDim(2) xDim(1) xDim(1)],[yDim(1) yDim(1) yDim(2) yDim(2) yDim(1)],'color',custom.Color(2,:))
line(hh(2),[xDim(1) xDim(1)],[0 1e6],'color',custom.Color(2,:))
line(hh(2),[xDim(2) xDim(2)],[0 1e6],'color',custom.Color(2,:))
line(hh(3),[yDim(1) yDim(1)],[0 1e6],'color',custom.Color(2,:))
line(hh(3),[yDim(2) yDim(2)],[0 1e6],'color',custom.Color(2,:))
% colorbar off
toc