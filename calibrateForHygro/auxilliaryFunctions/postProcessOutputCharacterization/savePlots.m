function savePlots(ctrl)
% Export all figures
if ctrl.plotFlag
    mkdir(horzcat(cd,ctrl.fileSep,'Plots'))
    FigIdx = findall(0, 'type', 'figure');
    for mm = 1:length(FigIdx)
        figure(FigIdx(mm))
        foo = get(gcf,'Name');
        foo = strrep(foo, '.', ',');
        print(foo,'-dpng','-r600')    
        movefile(horzcat(foo,'.png'),horzcat(cd,ctrl.fileSep,'Plots',ctrl.fileSep,foo,'.png'))
    end
end