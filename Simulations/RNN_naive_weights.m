
clc
clear all
close all

%%

load('./files/avwin.mat')

figure(1)
colormap('parula')
imagesc(avwin)
box off
cb = colorbar;
set(cb, 'tickdir', 'out', 'ytick', [0:0.02:0.1])
caxis([0, 0.1])
set(gca, 'ytick', [1:8], 'xtick', [0 3 6 7]+0.5, 'yticklabel', [], 'xticklabel', [], 'tickdir', 'out','fontsize', 18, 'linewidth', 1.5)
set(gca,'TickLength',[0.02, 0.01])

FigW = 4.5 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/w0in.pdf')

%%

load('./files/avwrec.mat')

figure(2)
colormap('parula')
imagesc(avwrec)
box off
cb = colorbar;
set(cb, 'tickdir', 'out', 'ytick', [0:0.02:0.06])
caxis([0, 0.06])
set(gca, 'ytick', [1:8], 'xtick', [1:8], 'yticklabel', [], 'xticklabel', [], 'tickdir', 'out','fontsize', 18, 'linewidth', 1.5)
set(gca,'TickLength',[0.02, 0.01])

FigW = 4.7 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/w0rec.pdf')

%%

load('./files/avwout.mat')

figure(3)
colormap('parula')
imagesc(avwout(1:4))
box off
caxis([0, 0.03])
set(gca, 'ytick', [1:8], 'xtick', [1:8], 'yticklabel', [], 'xticklabel', [], 'tickdir', 'out','fontsize', 18, 'linewidth', 1.5)
set(gca,'TickLength',[0.02, 0.01])

FigH = 4 ;
FigW = 0.5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/w0out.pdf')

figure(4)
colormap('parula')
cb = colorbar;
set(cb, 'tickdir', 'out','ytick', [0:0.01:0.03])
set(gca, 'ytick', [1:8], 'xtick', [1:8], 'yticklabel', [], 'xticklabel', [], 'tickdir', 'out','fontsize', 18, 'linewidth', 1.5)
caxis([0, 0.03])
set(gca,'TickLength',[0.02, 0.01])

FigH = 4 ;
FigW = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/cb.pdf')