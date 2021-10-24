clc
clear all
close all

%%

% reward probabilities
prob_mdprl(:,:,1)           = [0.92, 0.75, 0.43 ;
                               0.50, 0.50, 0.50 ;
                               0.57, 0.25, 0.08];
prob_mdprl(:,:,2)           = [0.16, 0.75, 0.98 ;
                               0.50, 0.50, 0.50 ;
                               0.02, 0.25, 0.84];
prob_mdprl(:,:,3)           = [0.92, 0.75, 0.43 ;
                               0.50, 0.50, 0.50 ;
                               0.57, 0.25, 0.08];

% assigning a unique index to RL agents following value of features/conjucntions/objects
for d = 1:3
    index_shp(:,:,d)     = [0, 1, 2; 0, 1, 2; 0, 1, 2] ;
    index_pttrn(:,:,d)   = [0, 0, 0; 1, 1, 1; 2, 2, 2] ;
    index_clr(:,:,d)     = [1, 1, 1; 1, 1, 1; 1, 1, 1]*d ;
    index_shppttrn(:,:,d)= index_shp(:,:,d)*3   + index_pttrn(:,:,d) ;
    index_pttrnclr(:,:,d)= index_pttrn(:,:,d)*3 + (index_clr(:,:,d)-1) ;
    index_shpclr(:,:,d)  = index_shp(:,:,d)*3   + (index_clr(:,:,d)-1) ;
end
index_shp = index_shp + 1 ;
index_pttrn = index_pttrn + 4 ;
index_clr = index_clr + 6 ;
index_shppttrn = index_shppttrn + 10 ;
index_pttrnclr = index_pttrnclr + 19 ;
index_shpclr = index_shpclr + 28 ;

%%

meanv   = inline('mean(x(:))','x') ;
vec     = inline('x(:)','x') ;
nrep    = 10 ;
lrm     = 0.05 ;
nTrials     = 12*27 ;
nTrialsL    = 27 ;
dc          = 0.005 ;

DA = zeros(nTrials,nrep) ;          % reward assignment
pop = zeros(nTrials,nrep) ;         % input object
v = 0.5*ones(nTrials,63,nrep) ;     % value traces
lr = [lrm*ones(1,9) lrm*ones(1,27) lrm*ones(1,27)] ;    % learning rates
for cntRep = 1:nrep
    for cnt_sd = 1:nTrials/27
        popidx = randperm(27) ;
        DA((cnt_sd-1)*27+1:cnt_sd*27,cntRep) = prob_mdprl(popidx)>rand(1,27) ;
        pop((cnt_sd-1)*27+1:cnt_sd*27,cntRep) = popidx ;
    end
    for cntTrials = 1:nTrials
       idx = [index_shp(pop(cntTrials,cntRep)) index_pttrn(pop(cntTrials,cntRep)) ...
              index_clr(pop(cntTrials,cntRep)) index_shppttrn(pop(cntTrials,cntRep)) ...
              index_pttrnclr(pop(cntTrials,cntRep)) index_shpclr(pop(cntTrials,cntRep)) ...
              pop(cntTrials,cntRep)+36] ;
       idxn = find(~ismember(1:63, idx)) ;
       v(cntTrials+1,idx,cntRep)  = v(cntTrials,idx,cntRep) + lr(idx).*(DA(cntTrials,cntRep)-v(cntTrials,idx,cntRep)) ;
       v(cntTrials+1,idxn,cntRep) = v(cntTrials,idxn,cntRep) + dc*(0.5-v(cntTrials,idxn,cntRep)) ;

       vFt(cntTrials,:,cntRep)    = v(cntTrials,index_pttrn,cntRep) ;
       vCnj(cntTrials,:,cntRep)   = v(cntTrials,index_shpclr,cntRep) ;
       vFtCnj(cntTrials,:,cntRep) = 0.5*(v(cntTrials,index_pttrn,cntRep)+v(cntTrials,index_shpclr,cntRep)) ;
       vObj(cntTrials,:,cntRep)   = v(cntTrials,63-27+1:63,cntRep) ;

       er(cntTrials,cntRep,1) = meanv((prob_mdprl-reshape(vFt(cntTrials,:,cntRep),3,3,3)).^2) ;
       er(cntTrials,cntRep,2) = meanv((prob_mdprl-reshape(vCnj(cntTrials,:,cntRep),3,3,3)).^2) ;
       er(cntTrials,cntRep,3) = meanv((prob_mdprl-reshape(vFtCnj(cntTrials,:,cntRep),3,3,3)).^2) ;
       er(cntTrials,cntRep,4) = meanv((prob_mdprl-reshape(vObj(cntTrials,:,cntRep),3,3,3)).^2) ;
    end
end

%%

close all
clrmat = colormap('lines(3)') ;
clrmat = clrmat(:,[1 3 2]) ;

uX                  = [1:324] ;
mu1                 = nanmean(squeeze(er(:,:,1)),2)' ;
sd1                 = nanstd(squeeze(er(:,:,1)),[],2)'./sqrt(nrep) ;

mu2                 = nanmean(squeeze(er(:,:,2)),2)' ;
sd2                 = nanstd(squeeze(er(:,:,2)),[],2)'./sqrt(nrep) ;

mu3                 = nanmean(squeeze(er(:,:,3)),2)' ;
sd3                 = nanstd(squeeze(er(:,:,3)),[],2)'./sqrt(nrep) ;

mu4                 = nanmean(squeeze(er(:,:,4)),2)' ;
sd4                 = nanstd(squeeze(er(:,:,4)),[],2)'./sqrt(nrep) ;

% prepare for patch
x1                  = [uX fliplr(uX)];
y1                  = [mu1+sd1 fliplr(mu1-sd1)];

x2                  = [uX fliplr(uX)];
y2                  = [mu2+sd2 fliplr(mu2-sd2)];

x3                  = [uX fliplr(uX)];
y3                  = [mu3+sd3 fliplr(mu3-sd3)];

x4                  = [uX fliplr(uX)];
y4                  = [mu4+sd4 fliplr(mu4-sd4)];

figure(1)
hold on
hpatch    = patch(x1,y1,'k'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',clrmat(:,1)); 
hline     = plot(uX,mu1,'-','Color',clrmat(:,1)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,1)); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x2,y2,'k'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor','m'); 
hline     = plot(uX,mu2,'-','Color','m'); 
set(hline,'LineWidth',2); 
set(hline,'Color','m'); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x3,y3,'k'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',clrmat(:,3)); 
hline     = plot(uX,mu3,'-','Color',clrmat(:,3)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,3)); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x4,y4,'k'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',clrmat(:,2)); 
hline     = plot(uX,mu4,'-','Color',clrmat(:,2)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,2)); 
box off
alpha(hpatch,0.3);

set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'XTick',[1 100:100:300],...
        'ytick',[0.03:0.01:0.07])
set(gca,'TickDir','out')
ylabel('MSE')
xlabel('trial (within a session)')
axis([1 nTrials 0.03 0.07])

cd ./figures
FigW = 6;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','MSE.pdf')
cd ../
