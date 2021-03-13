
%##########################################################################
%                           Ömer KARAGÖZ 2017                             %
%                         okaragoz@ogu.edu.tr                             %
%                                                                         %
%               DEFORMED & UNDERFORMED SYSTEM PLOTTER                     %
%                      WITH COLORED STRESSESS                             %
%This code was written on 01.05.2017 within the scope of METU CE583 course%
%                                                                         %
%INSTRUCTOR: Dr. Ozgur KURC                                               %
%MATLAB R2015b - 64bit is used                                            %
%#########################################################################%
close all
figure
axis off
axis equal
normalizedu = (u/(max(abs(u))))*0.3;
for i=1:numElem
    tempConn = conn(i,:);
    tempElemCoord = nodeCoordElemList(i).nodeCoordElem;
    defElemCoord = [];
    for j=1:4
        nodeID = tempConn(j);
        defElemCoord = [defElemCoord; tempElemCoord(j,1)+...
         normalizedu(2*nodeID-1) tempElemCoord(j,2)+normalizedu(2*nodeID)];
    end
    nodeCoordElemList(i).nodeCoordElemDef = defElemCoord;
end
sigmatype = 1;          % 1-> SigmaX, 2-> SigmaY, 3-> TaoXY
ssigma = [sigmaElemList.sigma];
maxsigma = max(max(ssigma(sigmatype,:)));
minsigma = min(min(ssigma(sigmatype,:)));
rangesigma = (maxsigma-minsigma);
dsigma = rangesigma/10;
colordata = jet(ceil(rangesigma));
for i=1:numElem
    cursigmamax = max(sigmaElemList(i).sigma(sigmatype,:));
    cursigmamin = min(sigmaElemList(i).sigma(sigmatype,:));
    if (abs(cursigmamin) <= abs(cursigmamax+1E-5))
        cursigma = cursigmamax;
    else
        cursigma = cursigmamin;
    end
    c = colordata(ceil(abs(minsigma-cursigma+1E-5)),:);
    tempElemCoord = nodeCoordElemList(i).nodeCoordElemDef;
    patch([tempElemCoord(1,1) tempElemCoord(2,1) tempElemCoord(3,1)...
        tempElemCoord(4,1)],[tempElemCoord(1,2) tempElemCoord(2,2)...
        tempElemCoord(3,2) tempElemCoord(4,2)],...
        [c(1) c(2) c(3)],'EdgeColor','none');
end
for i=1:numElem
    tempElemCoord = nodeCoordElemList(i).nodeCoordElem;
    patch([tempElemCoord(1,1) tempElemCoord(2,1) tempElemCoord(3,1)...
        tempElemCoord(4,1)],[tempElemCoord(1,2) tempElemCoord(2,2)...
        tempElemCoord(3,2) tempElemCoord(4,2)],...
        'green','FaceAlpha',0,'LineWidth',1);
end
hcb = colorbar('eastoutside');
colormap(jet(ceil(rangesigma)))
set(hcb,'Ticks',[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
tickLabels = {};
for i=0:9
    tickLabels = [tickLabels, num2str(sprintf('%10.2f',minsigma+i*dsigma))];
end
tickLabels = [tickLabels, sprintf('%10.2f',maxsigma)];
set(hcb,'YTickLabel',tickLabels)
saveas(gcf,'stress.bmp');
