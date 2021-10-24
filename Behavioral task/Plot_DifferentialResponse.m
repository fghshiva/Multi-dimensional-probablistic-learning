
clc
clear
close all
% rng('shuffle')
randstate = clock ;

%%

load('./files/RPL2Analysisv3_5_FeatureBased') ;
obj         = load('./files/RPL2Analysisv3_5_ObjectBased') ;
conj        = load('./files/RPL2Analysisv3_5_ConjunctionBased') ;

test        = str2func('signrank') ;
idxSubject  = [1:length(subjects)] ;

%%

for cnt_sbj = 1:length(subjects)
    inputname   = ['./PRLexp/inputs/input_', subjects{cnt_sbj} , '.mat'] ;
    resultsname = ['./PRLexp/SubjectData/PRL_', subjects{cnt_sbj} , '.mat'] ;

    load(inputname)
    load(resultsname)

    ntrialPerf          = [33:length(results.reward)] ;%[33:432] ;
    perfTH              = 0.5 + 2*sqrt(.5*.5/length(ntrialPerf)) ;
    [~, idxMax]         = max(expr.prob{1}(input.inputTarget)) ;
    choiceRew{cnt_sbj}  = results.choice' == idxMax' ;
    perf(cnt_sbj)       = nanmean(choiceRew{cnt_sbj}(33:432)) ;
    
    expr.shapeMap       = expr.targetShape ;
    expr.colorMap       = expr.targetColor ;
    expr.patternMap     = expr.targetPattern ;
    
    input.inputTarget   = input.inputTarget(:, ntrialPerf) ;
    expr.Ntrials        = length(ntrialPerf) ;
    results.reward      = results.reward(ntrialPerf) ;
    results.choice      = results.choice(ntrialPerf) ;

    idx_rew             = find(results.reward(2:end-1)==1)+1 ;
    idx_unr             = find(results.reward(2:end-1)==0)+1 ;
    targCh              = input.inputTarget(results.choice'+2*(0:(expr.Ntrials-1))) ;
    
    pFtTrials{cnt_sbj}      = nan*ones(1, length(ntrialPerf)) ;
    pFt0Trials{cnt_sbj}     = nan*ones(1, length(ntrialPerf)) ;
    
    pObTrials{cnt_sbj}      = nan*ones(1, length(ntrialPerf)) ;
    pOb0Trials{cnt_sbj}     = nan*ones(1, length(ntrialPerf)) ;
    
    pFtinfTrials{cnt_sbj}   = nan*ones(1, length(ntrialPerf)) ;
    pFtinf0Trials{cnt_sbj}  = nan*ones(1, length(ntrialPerf)) ;
    
    pFtninfTrials{cnt_sbj}  = nan*ones(1, length(ntrialPerf)) ;
    pFtninf0Trials{cnt_sbj} = nan*ones(1, length(ntrialPerf)) ;
    
    pCjTrials{cnt_sbj}      = nan*ones(1, length(ntrialPerf)) ;
    pCj0Trials{cnt_sbj}     = nan*ones(1, length(ntrialPerf)) ;
    
    pCjnTrials{cnt_sbj}     = nan*ones(1, length(ntrialPerf)) ;
    pCjn0Trials{cnt_sbj}    = nan*ones(1, length(ntrialPerf)) ;
    
    for cnt_O = expr.playcombinations
        shapeO          = find(ismember(1:3, expr.shapeMap(cnt_O))) ;
        colorO          = find(ismember(1:3, expr.colorMap(cnt_O))) ;
        patternO        = find(ismember(1:3, expr.patternMap(cnt_O))) ;

        if expr.flaginf==1
            infO        = shapeO ;
            expr.inf    = expr.shapeMap ;
            cnt_uinf    = find(expr.patternMap==patternO | expr.colorMap==colorO) ;
            cnt_Tcj     = find(expr.patternMap==patternO & expr.colorMap==colorO) ;
            cnt_ucj     = find((expr.patternMap==patternO & expr.shapeMap==shapeO) | (expr.colorMap==colorO & expr.shapeMap==shapeO)) ;
        elseif expr.flaginf==2
            infO        = patternO ;
            expr.inf    = expr.patternMap ;
            cnt_uinf    = find(expr.shapeMap==shapeO | expr.colorMap==colorO) ;
            cnt_Tcj     = find(expr.shapeMap==shapeO & expr.colorMap==colorO) ;
            cnt_ucj     = find((expr.shapeMap==shapeO & expr.patternMap==patternO) | (expr.colorMap==colorO & expr.patternMap==patternO)) ;
        end
        
        cnt_Tall        = find(expr.shapeMap==shapeO | expr.colorMap==colorO | expr.patternMap==patternO ) ;
        cnt_Tall        = cnt_Tall(~ismember(cnt_Tall, cnt_O)) ;
        cnt_Tinf        = find(expr.inf==infO) ;
        cnt_Tinf        = cnt_Tinf(~ismember(cnt_Tinf, cnt_O)) ;
        cnt_Tninf       = cnt_Tall(~ismember(cnt_Tall, cnt_Tinf)) ;
        cnt_Tinf        = cnt_Tinf(~ismember(cnt_Tinf, cnt_uinf)) ;
        cnt_Tcj         = cnt_Tcj(~ismember(cnt_Tcj, cnt_O)) ;
        cnt_ucj         = cnt_ucj(~ismember(cnt_ucj, cnt_O)) ;
        cnt_Tninf       = cnt_Tninf(~ismember(cnt_Tninf, cnt_Tcj)) ;
        cnt_Trest       = find(expr.shapeMap~=shapeO & expr.colorMap~=colorO & expr.patternMap~=patternO ) ;
        
        % 1) inf feature, 2) noninf feature, 3) inf conj, 4) object
        idxCell{1}      = cnt_Tinf ;
        idxCell{2}      = cnt_Tninf ;
        idxCell{3}      = cnt_Tcj ;
        idxCell{4}      = cnt_ucj ;
        idxCell{5}      = cnt_O ;
        
        % find trials with reward on object in trial (i) and no object on (i+1)
        idx_rewO        = idx_rew(targCh(idx_rew)==cnt_O & sum(ismember(input.inputTarget(:, idx_rew+1),  cnt_O))==0 & sum(ismember(input.inputTarget(:, idx_rew+1),  cnt_Trest))==1) ;
        idx_unrO        = idx_unr(targCh(idx_unr)==cnt_O & sum(ismember(input.inputTarget(:, idx_unr+1),  cnt_O))==0 & sum(ismember(input.inputTarget(:, idx_unr+1),  cnt_Trest))==1) ;
        
        % find trials with no target and no similar features for both options on trial (i+1):
        idx_rewOI       = idx_rewO(sum(ismember(input.inputTarget(:, idx_rewO+1),  cnt_Tinf))==1) ;
        idx_unrOI       = idx_unrO(sum(ismember(input.inputTarget(:, idx_unrO+1),  cnt_Tinf))==1) ;
        
        idx_rewON       = idx_rewO(sum(ismember(input.inputTarget(:, idx_rewO+1),  cnt_Tninf))==1) ;
        idx_unrON       = idx_unrO(sum(ismember(input.inputTarget(:, idx_unrO+1),  cnt_Tninf))==1) ;
        
        idx_rewOT       = idx_rewO(sum(ismember(input.inputTarget(:, idx_rewO+1),  cnt_Tall))==1) ;
        idx_unrOT       = idx_unrO(sum(ismember(input.inputTarget(:, idx_unrO+1),  cnt_Tall))==1) ;
        
        idx_rewOC       = idx_rewO(sum(ismember(input.inputTarget(:, idx_rewO+1),  cnt_Tcj))==1) ;
        idx_unrOC       = idx_unrO(sum(ismember(input.inputTarget(:, idx_unrO+1),  cnt_Tcj))==1) ;
        
        idx_rewIC       = idx_rewO(sum(ismember(input.inputTarget(:, idx_rewO+1),  cnt_ucj))==1) ;
        idx_unrIC       = idx_unrO(sum(ismember(input.inputTarget(:, idx_unrO+1),  cnt_ucj))==1) ;
        
        % find trials with reward on object in trial (i) and object on (i+1)
        idx_rewOb       = idx_rew(targCh(idx_rew)==cnt_O & sum(ismember(input.inputTarget(:, idx_rew+1),  cnt_O))==1) ;
        idx_unrOb       = idx_unr(targCh(idx_unr)==cnt_O & sum(ismember(input.inputTarget(:, idx_unr+1),  cnt_O))==1) ;
        
        % find trials with target and no similar features for both options on trial (i+1):
        idx_rewOb       = idx_rewOb(sum(ismember(input.inputTarget(:, idx_rewOb+1),  cnt_Tall))==0) ;
        idx_unrOb       = idx_unrOb(sum(ismember(input.inputTarget(:, idx_unrOb+1),  cnt_Tall))==0) ;
        
        % save through time
        pFtTrials{cnt_sbj}(idx_rewOT+1)         = ismember(targCh(idx_rewOT+1), cnt_Tall) ;
        pFt0Trials{cnt_sbj}(idx_unrOT+1)        = ismember(targCh(idx_unrOT+1), cnt_Tall) ;
        
        pObTrials{cnt_sbj}(idx_rewOb+1)         = ismember(targCh(idx_rewOb+1), cnt_O) ;
        pOb0Trials{cnt_sbj}(idx_unrOb+1)        = ismember(targCh(idx_unrOb+1), cnt_O) ;
        
        pFtinfTrials{cnt_sbj}(idx_rewOI+1)      = ismember(targCh(idx_rewOI+1), cnt_Tinf) ;
        pFtinf0Trials{cnt_sbj}(idx_unrOI+1)     = ismember(targCh(idx_unrOI+1), cnt_Tinf) ;
        
        pFtninfTrials{cnt_sbj}(idx_rewON+1)     = ismember(targCh(idx_rewON+1), cnt_Tninf) ;
        pFtninf0Trials{cnt_sbj}(idx_unrON+1)    = ismember(targCh(idx_unrON+1), cnt_Tninf) ;
        
        pCjTrials{cnt_sbj}(idx_rewOC+1)         = ismember(targCh(idx_rewOC+1), cnt_Tcj) ;
        pCj0Trials{cnt_sbj}(idx_unrOC+1)        = ismember(targCh(idx_unrOC+1), cnt_Tcj) ;
        
        pCjnTrials{cnt_sbj}(idx_rewIC+1)        = ismember(targCh(idx_rewIC+1), cnt_ucj) ;
        pCjn0Trials{cnt_sbj}(idx_unrIC+1)       = ismember(targCh(idx_unrIC+1), cnt_ucj) ;
    end
    pFtAll(cnt_sbj)     = nanmean(pFtTrials{cnt_sbj})    - nanmean(pFt0Trials{cnt_sbj}) ;
    pObAll(cnt_sbj)     = nanmean(pObTrials{cnt_sbj})    - nanmean(pOb0Trials{cnt_sbj}) ;
    pFtinfAll(cnt_sbj)  = nanmean(pFtinfTrials{cnt_sbj}) - nanmean(pFtinf0Trials{cnt_sbj}) ;
    pFtninfAll(cnt_sbj) = nanmean(pFtninfTrials{cnt_sbj})- nanmean(pFtninf0Trials{cnt_sbj}) ;
    pCjinfAll(cnt_sbj)  = nanmean(pCjTrials{cnt_sbj})    - nanmean(pCj0Trials{cnt_sbj}) ;
    pCjninfAll(cnt_sbj) = nanmean(pCjnTrials{cnt_sbj})   - nanmean(pCjn0Trials{cnt_sbj}) ;
    
    pFt(cnt_sbj,1)     = nanmean(pFtTrials{cnt_sbj}) ;
    pOb(cnt_sbj,1)     = nanmean(pObTrials{cnt_sbj}) ;
    pFtinf(cnt_sbj,1)  = nanmean(pFtinfTrials{cnt_sbj}) ;
    pFtninf(cnt_sbj,1) = nanmean(pFtninfTrials{cnt_sbj}) ;
    pCjinf(cnt_sbj,1)  = nanmean(pCjTrials{cnt_sbj}) ;
    pCjninf(cnt_sbj,1) = nanmean(pCjnTrials{cnt_sbj}) ;
    
    pFt(cnt_sbj,2)     = nanmean(pFt0Trials{cnt_sbj}) ;
    pOb(cnt_sbj,2)     = nanmean(pOb0Trials{cnt_sbj}) ;
    pFtinf(cnt_sbj,2)  = nanmean(pFtinf0Trials{cnt_sbj}) ;
    pFtninf(cnt_sbj,2) = nanmean(pFtninf0Trials{cnt_sbj}) ;
    pCjinf(cnt_sbj,2)  = nanmean(pCj0Trials{cnt_sbj}) ;
    pCjninf(cnt_sbj,2) = nanmean(pCjn0Trials{cnt_sbj}) ;
end

%%

idxperf             = (perf>perfTH) ;
idxperf(find(pObAll==0)) = 0 ;
idxperf             = find(idxperf) ;

%%

for cntDfit = 1:3
    RL2conj_couple{cntDfit}      = cat(1,conj.mlparRL2conj_couple{cntDfit,[idxSubject]}) ;
    RL2conj_uncouple{cntDfit}    = cat(1,conj.mlparRL2conj_uncouple{cntDfit,[idxSubject]}) ;
    RL2conj_decay{cntDfit}       = cat(1,conj.mlparRL2conj_decay{cntDfit,[idxSubject]}) ;
end

RL2ft_couple        = cat(1,mlparRL2_couple{[idxSubject]}) ;
RL2ft_uncouple      = cat(1,mlparRL2_uncouple{[idxSubject]}) ;
RL2ft_decay         = cat(1,mlparRL2_decay{[idxSubject]}) ;

RL2obj_couple       = cat(1,obj.mlparRL2obj_couple{[idxSubject]}) ;
RL2obj_uncouple     = cat(1,obj.mlparRL2obj_uncouple{[idxSubject]}) ;
RL2obj_decay        = cat(1,obj.mlparRL2obj_decay{[idxSubject]}) ;

N                   = [6 6 7 4 4 5 7 7 8] ;
Nt                  = 432 ;

LL                  = ([RL2ft_couple(:,100)          RL2ft_uncouple(:,100)            RL2ft_decay(:,100) ...
                        RL2obj_couple(:,100)         RL2obj_uncouple(:,100)           RL2obj_decay(:,100) ...
                        RL2conj_couple{1}(:,100)     RL2conj_uncouple{1}(:,100)       RL2conj_decay{1}(:,100) ...
                        RL2conj_couple{2}(:,100)     RL2conj_uncouple{2}(:,100)       RL2conj_decay{2}(:,100) ...
                        RL2conj_couple{3}(:,100)     RL2conj_uncouple{3}(:,100)       RL2conj_decay{3}(:,100)]) ;
                    
AIC                 = ([2*N(1)+2*RL2ft_couple(:,100)         2*N(2)+2*RL2ft_uncouple(:,100)           2*N(3)+2*RL2ft_decay(:,100) ...
                        2*N(4)+2*RL2obj_couple(:,100)        2*N(5)+2*RL2obj_uncouple(:,100)          2*N(6)+2*RL2obj_decay(:,100) ...
                        2*N(7)+2*RL2conj_couple{1}(:,100)    2*N(8)+2*RL2conj_uncouple{1}(:,100)      2*N(9)+2*RL2conj_decay{1}(:,100) ...
                        2*N(7)+2*RL2conj_couple{2}(:,100)    2*N(8)+2*RL2conj_uncouple{2}(:,100)      2*N(9)+2*RL2conj_decay{2}(:,100) ...
                        2*N(7)+2*RL2conj_couple{3}(:,100)    2*N(8)+2*RL2conj_uncouple{3}(:,100)      2*N(9)+2*RL2conj_decay{3}(:,100)]) ;

BIC                 = ([N(1)*log(Nt)+2*RL2ft_couple(:,100)         N(2)*log(Nt)+2*RL2ft_uncouple(:,100)             N(3)*log(Nt)+2*RL2ft_decay(:,100) ...
                        N(4)*log(Nt)+2*RL2obj_couple(:,100)        N(5)*log(Nt)+2*RL2obj_uncouple(:,100)          N(6)*log(Nt)+2*RL2obj_decay(:,100) ...
                        N(7)*log(Nt)+2*RL2conj_couple{1}(:,100)    N(8)*log(Nt)+2*RL2conj_uncouple{1}(:,100)      N(9)*log(Nt)+2*RL2conj_decay{1}(:,100) ...
                        N(7)*log(Nt)+2*RL2conj_couple{2}(:,100)    N(8)*log(Nt)+2*RL2conj_uncouple{2}(:,100)      N(9)*log(Nt)+2*RL2conj_decay{2}(:,100) ...
                        N(7)*log(Nt)+2*RL2conj_couple{3}(:,100)    N(8)*log(Nt)+2*RL2conj_uncouple{3}(:,100)      N(9)*log(Nt)+2*RL2conj_decay{3}(:,100)]) ;

[idxFtMd idxObjMd idxConjMd idxModel] = fIndex(LL) ;

idxFt       = idxFtMd ;
idxObj      = idxObjMd ;
idxConj     = idxConjMd{1} ;

idxFt       = idxperf(ismember(idxperf, idxFt)) ;
idxObj      = idxperf(ismember(idxperf, idxObj)) ;
idxConj     = idxperf(ismember(idxperf, idxConj)) ;

%%

funcM = @(x) nanmean(x) ;
range = [-0.6:0.1:0.6] ;
figure(19)
hold on
plot(pFtinf(idxperf,2),pFtinf(idxperf,1), 'db', 'linewidth', 2, 'markersize', 8)
plot(0:0.1:1, 0:0.1:1, '--k')
axis([0 1 0 1])
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',0:0.2:1)
box off
set(gca, 'tickdir', 'out')

axes('position', [0.65 0.25 0.2 0.2])
hold on
histogram(pFtinf(idxperf,1)-pFtinf(idxperf,2), range, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.7) ;
plot(funcM(pFtinf(idxperf,1)-pFtinf(idxperf,2))*ones(10,1), linspace(0,19.6, 10),'--b', 'linewidth', 1)
plot(0*ones(10,1), linspace(0,19.6, 10),'-', 'color', 0*[1 1 1], 'linewidth', 1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:10:20,'Xtick',-0.6:0.6:0.6)
box off
set(gca, 'tickdir', 'out')
axis([-0.6 0.6 0 20])

% cd ./figures
FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','pS_ft.pdf')
% cd ../

figure(20)
hold on
plot(pFtninf(idxperf,2),pFtninf(idxperf,1), 'db', 'linewidth', 2, 'markersize', 8)
plot(0:0.1:1, 0:0.1:1, '--k')
axis([0 1 0 1])
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',0:0.2:1)
box off
set(gca, 'tickdir', 'out')

axes('position', [0.65 0.25 0.2 0.2])
hold on
histogram(pFtninf(idxperf,1)-pFtninf(idxperf,2), range, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.7) ;
plot(funcM(pFtninf(idxperf,1)-pFtninf(idxperf,2))*ones(10,1), linspace(0,19.6, 10),'--b', 'linewidth', 1)
plot(0*ones(10,1), linspace(0,19.6, 10),'-', 'color', 0*[1 1 1], 'linewidth', 1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:10:20,'Xtick',-0.6:0.6:0.6)
box off
set(gca, 'tickdir', 'out')
axis([-0.6 0.6 0 20])

% cd ./figures
FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','pS_nft.pdf')
% cd ../

figure(21)
hold on
plot(pCjinf(idxperf,2),pCjinf(idxperf,1), 'db', 'linewidth', 2, 'markersize', 8)
plot(0:0.1:1, 0:0.1:1, '--k')
axis([0 1 0 1])
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',0:0.2:1)
box off
set(gca, 'tickdir', 'out')

axes('position', [0.65 0.25 0.2 0.2])
hold on
histogram(pCjinf(idxperf,1)-pCjinf(idxperf,2), range, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.7) ;
plot(funcM(pCjinf(idxperf,1)-pCjinf(idxperf,2))*ones(10,1), linspace(0,19.6, 10),'--b', 'linewidth', 1)
plot(0*ones(10,1), linspace(0,19.6, 10),'-', 'color', 0*[1 1 1], 'linewidth', 1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:10:20,'Xtick',-0.6:0.6:0.6)
box off
set(gca, 'tickdir', 'out')
axis([-0.6 0.6 0 20])

% cd ./figures
FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','pS_cj.pdf')
% cd ../

figure(22)
hold on
plot(pCjninf(idxperf,2),pCjninf(idxperf,1), 'db', 'linewidth', 2, 'markersize', 8)
plot(0:0.1:1, 0:0.1:1, '--k')
axis([0 1 0 1])
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:0.2:1,'Xtick',0:0.2:1)
box off
set(gca, 'tickdir', 'out')

axes('position', [0.65 0.25 0.2 0.2])
hold on
histogram(pCjninf(idxperf,1)-pCjninf(idxperf,2), range, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.7) ;
plot(funcM(pCjninf(idxperf,1)-pCjninf(idxperf,2))*ones(10,1), linspace(0,19.6, 10),'--b', 'linewidth', 1)
plot(0*ones(10,1), linspace(0,19.6, 10),'-', 'color', 0*[1 1 1], 'linewidth', 1)
set(gca,'FontName','Helvetica','FontSize',23,'FontWeight','normal','LineWidth',2,'yTick',0:10:20,'Xtick',-0.6:0.6:0.6)
box off
set(gca, 'tickdir', 'out')
axis([-0.6 0.6 0 20])

% cd ./figures
FigW = 6 ;
FigH = 5 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',[0,0,FigW,FigH],'units','centimeters');  
% print('-dpdf','-r500','pS_ncj.pdf')
% cd ../