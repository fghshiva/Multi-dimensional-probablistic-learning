
clc
clear 
close all

%%

load('./files/RPL2Analysisv3_5_FeatureBased') ;
obj              = load('./files/RPL2Analysisv3_5_ObjectBased') ;
conj             = load('./files/RPL2Analysisv3_5_ConjunctionBased') ;

test             = str2func('signrank') ;
cntDfit          = 1 ;
Nrep             = 10 ;

ntrialPerf       = [33:432] ;
perfTH           = 0.5 + 2*sqrt(.5*.5/length(ntrialPerf)) ;
idxSubject       = [1:length(subjects)] ;

probeTrialsAll   = expr.trialProbe ;
wSize            = 20 ; 

%%

RL2conj_couple      = cat(1,conj.mlparRL2conj_couple{cntDfit,idxSubject}) ;
RL2conj_uncouple    = cat(1,conj.mlparRL2conj_uncouple{cntDfit,idxSubject}) ;
RL2conj_decay       = cat(1,conj.mlparRL2conj_decay{cntDfit,idxSubject}) ;

RL2ft_couple        = cat(1,mlparRL2_couple{idxSubject}) ;
RL2ft_uncouple      = cat(1,mlparRL2_uncouple{idxSubject}) ;
RL2ft_decay         = cat(1,mlparRL2_decay{idxSubject}) ;

RL2obj_couple       = cat(1,obj.mlparRL2obj_couple{idxSubject}) ;
RL2obj_uncouple     = cat(1,obj.mlparRL2obj_uncouple{idxSubject}) ;
RL2obj_decay        = cat(1,obj.mlparRL2obj_decay{idxSubject}) ;

%%

for cnt_probe = 1:length(probeTrialsAll)
    pEstAll{cnt_probe}   = nan*ones(length(subjects),27) ;
    XAll{cnt_probe}      = [] ;
end

clear sesdata
for cnt_sbj = 1:length(subjects)
    inputname   = ['./PRLexp/inputs/input_', subjects{cnt_sbj} , '.mat'] ;
    resultsname = ['./PRLexp/SubjectData/PRL_', subjects{cnt_sbj} , '.mat'] ;
    
    load(inputname)
    load(resultsname)
    
    sesdata.flagUnr               = 1 ;
    sesdata.sig                   = 0.2 ;
    sesdata.input                 = input ;
    sesdata.expr                  = expr ;
    sesdata.results               = results ;
    
    rew{cnt_sbj}                  = results.reward ;
    [~, idxMax]                   = max(expr.prob{1}(input.inputTarget)) ;
    choiceRew{cnt_sbj}            = results.choice' == idxMax ;
    perfMean(cnt_sbj)             = nanmean(choiceRew{cnt_sbj}(ntrialPerf)) ;

    %% RL2 feature decay
    expr.shapeMap               = repmat([1 2 3 ;
                                   1 2 3 ;
                                   1 2 3 ], 1,1,3) ;

    expr.colorMap                = repmat([1 1 1 ;
                                    2 2 2 ;
                                    3 3 3], 1,1,3)+3 ;
                
    expr.patternMap(:,:,1)        = ones(3,3)+6 ;
    expr.patternMap(:,:,2)        = 2*ones(3,3)+6 ;
    expr.patternMap(:,:,3)        = 3*ones(3,3)+6 ;
    sesdata.expr                  = expr ;
    
    sesdata.flag_couple     = 0 ;
    sesdata.flag_updatesim  = 0 ;
    NparamBasic             = 5 ;
    if sesdata.flagUnr==1
        sesdata.Nalpha      = 2 ;
    else
        sesdata.Nalpha      = 1 ;
    end
    
    xpar                             = RL2ft_decay(cnt_sbj, 1:NparamBasic+sesdata.Nalpha) ;
    LLAll_RL2(cnt_sbj,:)             = fMLchoiceLL_RL2decay(xpar, sesdata) ;

    %% RL2 object decay
    sesdata.flag_couple = 0 ;
    NparamBasic         = 3 ;
    if sesdata.flagUnr==1
        sesdata.Nalpha  = 2 ;
    else
        sesdata.Nalpha  = 1 ;
    end
    xpar                            = RL2obj_decay(cnt_sbj, 1:NparamBasic+sesdata.Nalpha) ;
    LLAll_RL2obj(cnt_sbj,:)         = fMLchoiceLL_RL2objdecay(xpar, sesdata) ;
        
    %% RL2 conj decay
    expr.shapeMap               = repmat([1 2 3 ;
                                   1 2 3 ;
                                   1 2 3 ], 1,1,3) ;

    expr.colorMap                = repmat([1 1 1 ;
                                    2 2 2 ;
                                    3 3 3], 1,1,3) ;
                                
    expr.patternMap(:,:,1)        = ones(3,3) ;
    expr.patternMap(:,:,2)        = 2*ones(3,3) ;
    expr.patternMap(:,:,3)        = 3*ones(3,3) ;
    sesdata.expr                  = expr ;
    
    sesdata.cntD        = cntDfit ;
    sesdata.flag_couple = 0 ;
    NparamBasic         = 4 ;
    if sesdata.flagUnr==1
        sesdata.Nalpha = 4 ;
    else
        sesdata.Nalpha = 2 ;
    end
    xpar                            = RL2conj_decay(cnt_sbj, 1:NparamBasic+sesdata.Nalpha) ;
    LLAll_RL2conj(cnt_sbj,:)        = fMLchoiceLL_RL2conjdecay(xpar, sesdata) ;
    
    %% reported estimates
        
    for cnt_probe = 1:length(results.probEst)
        % read the estimated Odds ratio and convert to probability
        probEst         = results.probEst{cnt_probe} ;
        probEst         = probEst./(1+probEst) ;
        
        % regression coefficient
        X   = 4 ; 
        ll1 = [X      1  1/X;
               1/X^2  1  X^2;
               X      1  1/X]' ;

        X   = 3 ; 
        ll  = [X      1  1/X]' ;
        
        LLSh            = nan*ones(3,3,3) ;
        LLSh(1,:,:)     = ll1*ll(1) ;
        LLSh(2,:,:)     = ones(size(ll1))*ll(2) ;
        LLSh(3,:,:)     = ll1*ll(3) ;
        
        % informative feature
        llRegft = [1  1  1;
                   1  1  1;
                   1  1  1]' ;
        LLRegft          = nan*ones(3,3,3) ;
        LLRegft(1,:,:)   = llRegft*ll(1) ;
        LLRegft(2,:,:)   = llRegft*ll(2) ;
        LLRegft(3,:,:)   = llRegft*ll(3) ;
        Regft            = LLRegft./(1+LLRegft) ;
        
        % conjunction
        llRegcnj = [X      1  1/X;
                    1/X^2  1  X^2;
                    X      1  1/X]' ;
        LLRegcnj         = nan*ones(3,3,3) ;
        LLRegcnj(1,:,:)  = llRegcnj.^(2/3)*1 ;
        LLRegcnj(2,:,:)  = llRegcnj.^(2/3)*1 ;
        LLRegcnj(3,:,:)  = llRegcnj.^(2/3)*1 ;
        Regcnj           = LLRegcnj./(1+LLRegcnj) ;
        
         % objects
        LLRegobj         = LLSh ;
        Regobj           = LLRegobj./(1+LLRegobj) ;
        
        %%
        if (perfMean(cnt_sbj)>perfTH)
            pEstAll{cnt_probe}(cnt_sbj,:)      = nan*probEst(expr.playcombinations) ;
            XTEMP                              = [Regft(expr.playcombinations); Regcnj(expr.playcombinations); Regobj(expr.playcombinations); probEst(expr.playcombinations)]' ; 
            XTEMP                              = round(XTEMP./0.05)*0.05 ;
            XAll{cnt_probe}                    = [XAll{cnt_probe}; XTEMP] ;
        else
            pEstAll{cnt_probe}(cnt_sbj,:)      = nan*probEst(expr.playcombinations) ;
            XTEMP                              = [Regft(expr.playcombinations); Regcnj(expr.playcombinations); Regobj(expr.playcombinations); nan*probEst(expr.playcombinations)]' ; 
            XTEMP                              = round(XTEMP./0.05)*0.05 ;
            XAll{cnt_probe}                    = [XAll{cnt_probe}; nan*XTEMP] ;
        end
    end
end

%%

idxperf             = (perfMean>perfTH) ;
idxperf(46)         = 0 ;
idxperf(54)         = 0 ;
idxperf             = find(idxperf) ;

%%
                    
clear b1 b2 b3 b4 b5
for cnt_sbj = 1:length(subjects)
    for cnt_probe = 1:length(probeTrialsAll)
        idxS         = (cnt_sbj-1)*length(expr.playcombinations)+1:(cnt_sbj)*length(expr.playcombinations) ;
        X            = XAll{cnt_probe}(idxS,end) ;
        Y1           = XAll{cnt_probe}(idxS,1) ;
        Y1           = [Y1 ones(size(Y1,1),1)] ;
        Y2           = XAll{cnt_probe}(idxS,[1 2]) ;
        Y2           = [Y2 ones(size(Y2,1),1)] ;
        Y3           = XAll{cnt_probe}(idxS,3) ;
        Y3           = [Y3 ones(size(Y3,1),1)] ;

        b1(:,cnt_probe,cnt_sbj)     = regress(X,Y1) ;
        b2(:,cnt_probe,cnt_sbj)     = regress(X,Y2) ;
        b3(:,cnt_probe,cnt_sbj)     = regress(X,Y3) ;
        yCalc1       = Y1*b1(:,cnt_probe,cnt_sbj) ;
        yCalc2       = Y2*b2(:,cnt_probe,cnt_sbj) ;
        yCalc3       = Y3*b3(:,cnt_probe,cnt_sbj) ;
        y            = X ;
        RsqS(1, cnt_probe,cnt_sbj)   = 1 - nansum((y - yCalc1).^2)/nansum((y - nanmean(y)).^2) ;
        RsqS(2, cnt_probe,cnt_sbj)   = 1 - nansum((y - yCalc2).^2)/nansum((y - nanmean(y)).^2) ;
        RsqS(3, cnt_probe,cnt_sbj)   = 1 - nansum((y - yCalc3).^2)/nansum((y - nanmean(y)).^2) ;
        
    end
    fvaltmp1    = 1 ;
    fvaltmp2    = 1 ;
    for cnt_rep = 1:Nrep
        x0                                     = [10*rand(1), 10^4*rand(1)];
        lb                                     = [0, 0]; 
        ub                                     = [10, 10^4];
        [parsEtmp1, resnorm1]                  = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,b2(1,:,cnt_sbj),lb,ub);
        [parsEtmp2, resnorm2]                  = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,b2(2,:,cnt_sbj),lb,ub);
        if  resnorm1 < fvaltmp1
            parEs(1, cnt_sbj, :) = parsEtmp1 ;
            fvaltmp1             = resnorm1 ;
        end
        if  resnorm2 < fvaltmp2
            parEs(2, cnt_sbj, :) = parsEtmp2 ;
            fvaltmp2             = resnorm2 ;
        end
    end
end

for cnt_rep = 1:Nrep
    % fit an exponential to data
    x0                                     = [1*rand(1), 1000*rand(1)];
    lb                                     = [0, 0]; 
    ub                                     = [1, 1000];
    [parEtmp1, resnorm1]                   = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,nanmedian(b2(1,:,idxperf),3),lb,ub);
    [parEtmp2, resnorm2]                   = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,nanmedian(b2(2,:,idxperf),3),lb,ub);
    
    fvaltmp1    = 1 ;
    fvaltmp2    = 1 ;
    if  resnorm1 < fvaltmp1
        parE(1, :)           = parEtmp1 ;
        fvaltmp1             = resnorm1 ;
    end
    if  resnorm2 < fvaltmp2
        parE(2, :)           = parEtmp2 ;
        fvaltmp2             = resnorm2 ;
    end
end

%%

close all
clrmat = colormap('lines(3)') ;
clrmat = clrmat(:,[1 3 2]) ;

clear X1 X2 X3 X4
dcshift  =0.6 ;
LL1                 = cat(2, LLAll_RL2(idxperf,:)') ;
LL1                 = nanmean(reshape(LL1, [], 1, length(idxperf)),2) ;
LL1                 = reshape(LL1,[], length(idxperf)) ;

LL2                 = cat(2, LLAll_RL2obj(idxperf,:)') ;
LL2                 = nanmean(reshape(LL2, [], 1, length(idxperf)),2) ;
LL2                 = reshape(LL2,[], length(idxperf)) ;

LL3                 = cat(2, LLAll_RL2conj(idxperf,:)') ;
LL3                 = nanmean(reshape(LL3, [], 1, length(idxperf)),2) ;
LL3                 = reshape(LL3,[], length(idxperf)) ;

% filter data using a moving average box
uX                  = wSize:432 ;
for cntSubjects = 1:length(idxperf)
    X1(:,cntSubjects)     = movavgFilt(LL1(:,cntSubjects)', wSize, 'Left')' ;
end
X1                  = X1(wSize:end,:) ;
mu1                 = nanmean(2*X1+7/432,2)' ;
sd1                 = nanstd(2*X1')./sqrt(length(idxperf)) ;

for cntSubjects = 1:length(idxperf)
    X2(:,cntSubjects)     = movavgFilt(LL2(:,cntSubjects)', wSize, 'Left')' ;
end
X2                  = X2(wSize:end,:) ;
mu2                 = nanmean(2*X2+5/432,2)' ;
sd2                 = nanstd(2*X2')./sqrt(length(idxperf)) ;


for cntSubjects = 1:length(idxperf)
    X3(:,cntSubjects)     = movavgFilt(LL3(:,cntSubjects)', wSize, 'Left')' ;
end
X3                  = X3(wSize:end,:) ;
mu3                 = nanmean(2*X3+8/432,2)' ;
sd3                 = nanstd(2*X3')./sqrt(length(idxperf)) ;

for cntSubjects = 1:length(idxperf)
    X4(:,cntSubjects)     = dcshift + movavgFilt(2*LL3(:,cntSubjects)'-...
        2*LL1(:,cntSubjects)'+1/432, wSize, 'Left')' ;
end
X4                  = X4(wSize:end,:) ;
mu4                 = nanmean(X4,2)' ;
sd4                 = nanstd(2*X4')./sqrt(length(idxperf)) ;

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
hpatch    = patch(x1,y1,'b'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',clrmat(:,1)); 
hline     = plot(uX,mu1,'-','Color',clrmat(:,1)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,1)); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x2,y2,'r'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',clrmat(:,2)); 
hline     = plot(uX,mu2,'-','Color',clrmat(:,2)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,2)); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x3,y3,'m'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',clrmat(:,3)); 
hline     = plot(uX,mu3,'-','Color',clrmat(:,3)); 
set(hline,'LineWidth',2); 
set(hline,'Color',clrmat(:,3)); 
box off
alpha(hpatch,0.3);

hpatch = patch(x4,y4,'k'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',[0 .0 0]);
hold on;
hline = plot(uX,mu4,'k-'); 
set(hline,'LineWidth',2); 
set(hline,'Color','k'); 
box off
alpha(hpatch,0.3);

plot(1:432, dcshift*ones(1,432), '--k')
set(gca,'FontName','Helvetica','FontSize',25,'FontWeight','normal','LineWidth',2,'XTick',[1 100:100:450],...
        'ytick',[dcshift-0.2:0.2:dcshift 0.8:0.2:1.4], 'yticklabel',[-0.2:0.2:0 0.8:0.2:1.4])
set(gca,'TickDir','out')
ylabel('goodness-of-fit')
xlabel('trial (within a session)')
axis([1 450 dcshift-0.2 1.401])

% cd ./figures
FigW = 6;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','LL.pdf')
% cd ../

%%

rewAll          = cat(2, rew{idxperf}) ;
rewAll          = nanmean(reshape(rewAll, [], 1, length(idxperf)),2) ;
rewAll          = reshape(rewAll,[], length(idxperf)) ;

chrewAll        = cat(2, choiceRew{idxperf}) ;
chrewAll        = nanmean(reshape(chrewAll, [], 1, length(idxperf)),2) ;
chrewAll        = reshape(chrewAll,[], length(idxperf)) ;

uX              = wSize:432 ;
X1              = filter(ones(1,wSize)/wSize,1,rewAll) ;
X1              = X1(wSize:end,:) ;
mu1             = nanmean(X1,2)' ;
sd1             = std(X1')./sqrt(length(idxperf)) ;

X2              = filter(ones(1,wSize)/wSize,1,chrewAll) ;
X2              = X2(wSize:end,:) ;
mu2             = nanmean(X2,2)' ;
sd2             = std(X2')./sqrt(length(idxperf)) ;

x               = [uX fliplr(uX)];
y1              = [mu1+sd1 fliplr(mu1-sd1)];
y2              = [mu2+sd2 fliplr(mu2-sd2)];

figure(2)
hpatch = patch(x,y1,'k'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',0*[1 1 1]); 
hold on;
hline = plot(uX,mu1,'k-'); 
set(hline,'LineWidth',2); 
set(hline,'Color','k'); 
box off
alpha(hpatch,0.3);

hpatch = patch(x,y2,[0.5 0.5 0.5]); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor',[0.5 0.5 0.5]); 
hold on;
hline = plot(uX,mu2,'-', 'Color', [0.5 0.5 0.5]); 
set(hline,'LineWidth',2); 
set(hline,'Color',[0.5 0.5 0.5]); 
box off
alpha(hpatch,0.3);

plot(1:432, 0.5*ones(1,432), '--k')

set(gca,'FontName','Helvetica','FontSize',25,'FontWeight','normal','LineWidth',2,'XTick',[1 100:100:450],'ytick',0:0.1:1)
for cntProbe = probeTrialsAll([1:5])
    x = [0.13+(cntProbe/450)*0.77 0.13+(cntProbe/450)*0.77];
    y = [0.22 0.32];
    annotation('textarrow',x,y)
end
set(gca,'TickDir','out')
axis([1 450 0.42 .7])
xlabel('trial (within a session)')
ylabel('performance')

% cd ./figures
FigW = 6;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','perf.pdf')
% cd ../

%%

clrmat = colormap('lines(3)') ;

figure(4)
for cnt_p = 1:3
    hold on
    errorbar(probeTrialsAll,nanmean(RsqS(cnt_p,:,idxperf),3),nanstd(RsqS(cnt_p,:,idxperf),[],3)./sqrt(length(idxperf)),'d', 'color', clrmat(:,cnt_p),'LineWidth',2, 'markersize', 8)
    set(gca,'FontName','Helvetica','FontSize',25,'FontWeight','normal','LineWidth',2,'yTick',0:0.1:0.4,'Xtick',[1 100:100:450], 'tickdir', 'out')
    box off
    axis([0 450 0. 0.4])
    
    for cnt_rep = 1:nrep
        % fit an exponential to data
        x0                                     = [10*rand(1), 10^4*rand(1)];
        lb                                     = [0, 0]; 
        ub                                     = [01, 10^4];
        [parEtmp, resnorm]                     = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,mean(RsqS(cnt_p,:,idxperf),3),lb,ub);

        fvaltmp    = 1 ;
        if  resnorm < fvaltmp
            parEr(1, :)          = parEtmp ;
            fvaltmp              = resnorm ;
        end
    end

    for cnt_s = 1:length(idxperf)
        fvaltmp    = 1 ;
        for cnt_rep = 1:nrep
            x0                                     = [10*rand(1), 10^4*rand(1)];
            lb                                     = [0, 0]; 
            ub                                     = [10, 10^4];
            [parsEtmp, resnorm]                    = lsqcurvefit(@Ffitexp, x0,probeTrialsAll,RsqS(cnt_p,:,idxperf(cnt_s)),lb,ub);
            if  resnorm < fvaltmp
                parEsr(1, cnt_s, :)  = parsEtmp ;
                fvaltmp             = resnorm ;
            end
        end
    end

    Xfit                                    = 1:450 ;
    Yfit1                                   = Ffitexp(parEr(1,:),Xfit) ;

    for cnt_s = 1:length(idxperf)
        YfitS(cnt_s,:)                  = Ffitexp(parEsr(1, cnt_s, :),Xfit) ;
    end
    Yfit1M              = Ffitexp(parEr(1, :),Xfit) ;

    uX                  = Xfit ;
    mu1                 = Yfit1 ;
    sd1                 = std(YfitS)./sqrt(length(idxperf)) ;
    
    x1                  = [uX fliplr(uX)];
    y1                  = [mu1+sd1 fliplr(mu1-sd1)];
    
    hpatch    = patch(x1,y1,'k'); 
    set(hpatch,'EdgeColor','none'); 
    set(hpatch,'FaceColor', clrmat(:,cnt_p)); 
    hline     = plot(uX,mu1); 
    set(hline,'LineWidth',2);
    set(hline,'Color',clrmat(:,cnt_p)); 
    box off
    alpha(hpatch,0.3);
    ylabel('R^2')
end
xlabel('trial (within a session)')

% cd ./figures
FigW = 6;
FigH = 5;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r10000','R2.pdf')
% cd ../

%%

figure(7)
hold on
errorbar(probeTrialsAll,nanmedian(b2(1,:,idxperf),3), nanstd(b2(1,:,idxperf),[],3)./sqrt(length(idxperf)), 's', 'color', clrmat(:,1),'LineWidth',2, 'markersize', 8)
errorbar(probeTrialsAll,nanmedian(b2(2,:,idxperf),3), nanstd(b2(2,:,idxperf),[],3)./sqrt(length(idxperf)), 's', 'color', [1 0 1],'LineWidth',2, 'markersize', 8)

Xfit                                    = 1:450 ;

Yfit1                                   = Ffitexp(parE(1,:),Xfit) ;
Yfit2                                   = Ffitexp(parE(2,:),Xfit) ;

for cnt_sbj = 1:length(subjects)
    Yfit1S(cnt_sbj,:)                  = Ffitexp(parEs(1, cnt_sbj, :),Xfit) ;
    Yfit2S(cnt_sbj,:)                  = Ffitexp(parEs(2, cnt_sbj, :),Xfit) ;
end
Yfit1M              = Ffitexp(parE(1, :),Xfit) ;
Yfit2M              = Ffitexp(parE(2, :),Xfit) ;

uX                  = Xfit ;
mu1                 = Yfit1M ;
sd1                 = std(Yfit1S(idxperf,:))./sqrt(length(idxperf)) ;

mu2                 = Yfit2M ;
sd2                 = std(Yfit2S(idxperf,:))./sqrt(length(idxperf)) ;

x1                  = [uX fliplr(uX)];
y1                  = [mu1+sd1 fliplr(mu1-sd1)];

x2                  = [uX fliplr(uX)];
y2                  = [mu2+sd2 fliplr(mu2-sd2)];

hpatch    = patch(x1,y1,'k'); 
set(hpatch,'EdgeColor','none'); 
set(hpatch,'FaceColor', clrmat(:,1)); 
hline     = plot(uX,mu1); 
set(hline,'LineWidth',2);
set(hline,'Color',clrmat(:,1)); 
box off
alpha(hpatch,0.3);

hpatch    = patch(x2,y2,'k'); 
set(hpatch,'EdgeColor','none');
set(hpatch,'FaceColor',[1 0 1]); 
hline     = plot(uX,mu2); 
set(hline,'LineWidth',2); 
set(hline,'Color',[1 0 1]); 
box off
alpha(hpatch,0.3);

set(gca,'FontName','Helvetica','FontSize',25,'FontWeight','normal','LineWidth',2,'yTick',0:0.1:1,'Xtick',[1 100:100:450], 'tickdir', 'out')
box off
axis([0 450 0. 0.4])
xlabel('trial (within a session)')
ylabel('normalized weight')

% cd ./figures
FigW = 6;
FigH = 5;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r10000','weights.pdf')
% cd ../
