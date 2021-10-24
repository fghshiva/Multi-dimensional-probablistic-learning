
clc
clear all
close all

%%

% reward probabilities
prob(:,:,1)           = [0.92, 0.75, 0.43 ;
                         0.50, 0.50, 0.50 ;
                         0.57, 0.25, 0.08];
prob(:,:,2)           = [0.16, 0.75, 0.98 ;
                         0.50, 0.50, 0.50 ;
                         0.02, 0.25, 0.84];
prob(:,:,3)           = [0.92, 0.75, 0.43 ;
                         0.50, 0.50, 0.50 ;
                         0.57, 0.25, 0.08];
                     
% assigning a unique index to features/conjucntions/objects
shapeMap          = repmat([1 2 3 ;
                            1 2 3 ;
                            1 2 3 ], 1,1,3) ;

colorMap          = repmat([1 1 1 ;
                            2 2 2 ;
                            3 3 3], 1,1,3) ;

patternMap(:,:,1) = ones(3,3) ;
patternMap(:,:,2) = 2*ones(3,3) ;
patternMap(:,:,3) = 3*ones(3,3) ;

shapecolorMap   = (shapeMap-1)*3+colorMap ;
shapepatternMap = (patternMap-1)*3+shapeMap ;
colorpattenrMap = (colorMap-1)*3+patternMap ;

%%

prbshp  = mean(mean(prob,1),3) ;
prbclr  = mean(mean(prob,3),2) ;
prbptrn = mean(mean(prob,1),2) ;

prbshpclr  = mean(prob,3) ;
prbshpptrn = mean(prob,1) ;
prbclrptrn = mean(prob,2) ;

probFt   = 0.3*(prbshp(shapeMap)+prbclr(colorMap)+prbptrn(patternMap)) ;
probCnj1 = prbshpclr(shapecolorMap) ;
probCnj2 = prbshpptrn(shapepatternMap) ;
probCnj3 = prbclrptrn(colorpattenrMap) ;

probFtCnj1 = 0.5*(prbptrn(patternMap)+prbshpclr(shapecolorMap)) ;
probFtCnj2 = 0.5*(prbclr(colorMap)+prbshpptrn(shapepatternMap)) ;
probFtCnj3 = 0.5*(prbshp(shapeMap)+prbclrptrn(colorpattenrMap)) ;

genFt0    = corr(prob(:), probFt(:)) ;
genFtCnj0 = max([corr(prob(:), probFtCnj1(:)) corr(prob(:), probFtCnj2(:)) ...
                 corr(prob(:), probFtCnj3(:))]) ;
             
for cnt_Sh = 1:10000
   prob_Sh =  reshape(rand(1,27),3,3,3) ;
   prbshp = mean(mean(prob_Sh,1),3) ;
   prbclr = mean(mean(prob_Sh,3),2) ;
   prbptrn = mean(mean(prob_Sh,1),2) ;
   
   prbshpclr = mean(prob_Sh,3) ;
   prbshpptrn = mean(prob_Sh,1) ;
   prbclrptrn = mean(prob_Sh,2) ;
   
   probFt_Sh = 0.3*(prbshp(shapeMap)+prbclr(colorMap)+prbptrn(patternMap)) ;
   probCnj1_Sh = prbshpclr(shapecolorMap) ;
   probCnj2_Sh = prbshpptrn(shapepatternMap) ;
   probCnj3_Sh = prbclrptrn(colorpattenrMap) ;
   
   probFtCnj1_Sh = 0.5*(prbptrn(patternMap)+prbshpclr(shapecolorMap)) ;
   probFtCnj2_Sh = 0.5*(prbclr(colorMap)+prbshpptrn(shapepatternMap)) ;
   probFtCnj3_Sh = 0.5*(prbshp(shapeMap)+prbclrptrn(colorpattenrMap)) ;
   
   genFt(cnt_Sh)    = corr(prob_Sh(:), probFt_Sh(:)) ;
   genFtCnj(cnt_Sh,:) = [corr(prob_Sh(:), probFtCnj1_Sh(:)) corr(prob_Sh(:), probFtCnj2_Sh(:)) ...
                           corr(prob_Sh(:), probFtCnj3_Sh(:))] ;
end
genFt    = genFt(:) ;
genFtCnj = genFtCnj(:) ;

%%

close all
range       = [-0.5:0.05:1] ;
xb          = 0.181 ;
figure(1)
hold on
histogram(genFt, range, 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'FaceColor', 'c','Normalization','probability')
plot(median(genFt)*ones(10,1), linspace(0,xb, 10), '--', 'Color', 'c','LineWidth',1)
plot(mean(genFt)*ones(10,1), linspace(0,xb, 10), '-', 'Color', 'c','LineWidth',1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal',...
    'LineWidth',2,'yTick',0:0.06:0.2,'Xtick',[-1:0.5:1], 'tickdir', 'out')
axis([range(1) range(end) 0 xb])
xlabel('generalizability')
ylabel('probability')

x = 0.175+0.73*(0.5+[genFt0 genFt0])/1.5 ;
y = [0.26 0.16]+0.03 ;
annotation('textarrow',x,y,'LineWidth',2)

FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','./figures/genFt.pdf')

%%

figure(2)
hold on
histogram(genFtCnj, range, 'EdgeColor', 'none', 'FaceAlpha', 0.65, 'FaceColor', [0.4470 0.3250 0.6940],'Normalization','probability')
plot(median(genFtCnj)*ones(10,1), linspace(0,xb, 10), '--', 'Color', [0.4470 0.3250 0.6940],'LineWidth',1)
plot(mean(genFtCnj)*ones(10,1), linspace(0,xb, 10), '-', 'Color', [0.4470 0.3250 0.6940],'LineWidth',1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal',...
    'LineWidth',2,'yTick',0:0.06:0.2,'Xtick',[-1:0.5:1], 'tickdir', 'out')
axis([range(1) range(end) 0 xb])
xlabel('generalizability')
ylabel('probability')

x = 0.175+0.73*(0.5+[genFtCnj0 genFtCnj0])/1.5 ;
y = [0.26 0.16]+0.03 ;
annotation('textarrow',x,y,'LineWidth',2)

FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','./figures/genFtCnj.pdf')
