function loglikehood = fMLchoicefit_RL2conj(xpar, sesdata)
%
% DESCRIPTION: fits data to RL(2) model using ML method
% 
% INPUT: 
% sesdata structure which includes input, experiment and behavioral data
% 
% OUTPUT:
% fitted parametres
% 
% Version History:
% 0.1:  [2015-09-13]

loglikehood = 0 ;
NparamBasic = 3 ;

xpar(2:3)= abs(xpar(2:3)) ;

BiasL = xpar(1) ;
magF  = xpar(2) ;
magC  = xpar(3) ;

xpar([NparamBasic+1:NparamBasic+sesdata.Nalpha])=1./(1+exp(-(xpar([NparamBasic+1:NparamBasic+sesdata.Nalpha]))./sesdata.sig) ) ;
alpha_rewColor      = xpar([NparamBasic+1]) ;
alpha_rewShape      = xpar([NparamBasic+1]) ;
alpha_rewPattern    = xpar([NparamBasic+1]) ; 
alpha_rew           = xpar([NparamBasic+2]) ;
if sesdata.flagUnr==1
    alpha_unrColor      = xpar([NparamBasic+3]) ;
    alpha_unrShape      = xpar([NparamBasic+3]) ;
    alpha_unrPattern    = xpar([NparamBasic+3]) ;
    alpha_unr           = xpar([NparamBasic+4]) ;
else
    alpha_unrColor      = alpha_rewColor ;
    alpha_unrShape      = alpha_rewShape ;
    alpha_unrPattern    = alpha_rewPattern ;
    alpha_unr           = alpha_rew ;
end

shapeMap        = sesdata.expr.shapeMap ;
colorMap        = sesdata.expr.colorMap ;
patternMap      = sesdata.expr.patternMap ;
inputTarget     = sesdata.input.inputTarget ;
correcttrials   = sesdata.results.reward ;
choicetrials    = sesdata.results.choice ;
flag_couple     = sesdata.flag_couple ;
flag_updatesim  = sesdata.flag_updatesim ;
ntrials         = length(choicetrials) ;
inputRewards    = sesdata.input.inputReward ;
cntD            = sesdata.cntD ;

vf              = (0.5*ones(3,1)) ; 
vc              = (0.5*ones(9,1)) ; 

for cnt_trial=1:ntrials
    
    correct         = correcttrials(cnt_trial) ;
    choice          = choicetrials(cnt_trial) ; 
    correctunCh     = inputRewards(3-choice, cnt_trial) ;
    choiceunCh      = 3-choice ;
    
    idx_shape(2)    = shapeMap(inputTarget(2, cnt_trial)) ;
    idx_color(2)    = colorMap(inputTarget(2, cnt_trial)) ;
    idx_pattern(2)  = patternMap(inputTarget(2, cnt_trial)) ;
    idx_shape(1)    = shapeMap(inputTarget(1, cnt_trial)) ;
    idx_color(1)    = colorMap(inputTarget(1, cnt_trial)) ;
    idx_pattern(1)  = patternMap(inputTarget(1, cnt_trial)) ;
    
    if cntD==1
        inputConj(1, cnt_trial) = (idx_pattern(1)-1)*3 + idx_shape(1) ;
        inputConj(2, cnt_trial) = (idx_pattern(2)-1)*3 + idx_shape(2) ;
        vp(1,:,:)       = magF*vf(1) + magC*vc ;
        vp(2,:,:)       = magF*vf(2) + magC*vc ;
        vp(3,:,:)       = magF*vf(3) + magC*vc ;
        pChoiceR = 1./(1+exp(-( (vp(inputTarget(2, cnt_trial))-vp(inputTarget(1, cnt_trial))) + BiasL ) )) ;
    elseif cntD==2
        inputConj(1, cnt_trial) = (idx_pattern(1)-1)*3 + idx_color(1) ;
        inputConj(2, cnt_trial) = (idx_pattern(2)-1)*3 + idx_color(2) ;
        vp(:,1,:)       = magF*vf(1) + magC*vc ;
        vp(:,2,:)       = magF*vf(2) + magC*vc ;
        vp(:,3,:)       = magF*vf(3) + magC*vc ;
        pChoiceR = 1./(1+exp(-( (vp(inputTarget(2, cnt_trial))-vp(inputTarget(1, cnt_trial))) + BiasL ) )) ;
    elseif cntD==3
        inputConj(1, cnt_trial) = (idx_shape(1)-1)*3 + idx_color(1) ;
        inputConj(2, cnt_trial) = (idx_shape(2)-1)*3 + idx_color(2) ;
        vp(:,:,1)       = magF*vf(1) + magC*vc ;
        vp(:,:,2)       = magF*vf(2) + magC*vc ;
        vp(:,:,3)       = magF*vf(3) + magC*vc ;
        pChoiceR = 1./(1+exp(-( (vp(inputTarget(2, cnt_trial))-vp(inputTarget(1, cnt_trial))) + BiasL ) )) ;
    end
    pChoiceL = 1-pChoiceR ;
    if cnt_trial >= 1  
        if choice == 2 
            loglikehood = loglikehood - log(pChoiceR) ;
        else
            loglikehood = loglikehood - log(pChoiceL) ; 
        end                      
    end
    
    % conjunction
    if correct
        idxC = inputConj(choice, cnt_trial) ;
        idxW = inputConj(3-choice, cnt_trial) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0) ;
        vc = update(vc, idxC, idxW, alpha_rew) ;
    else
        idxC = inputConj(3-choice, cnt_trial) ;
        idxW = inputConj(choice, cnt_trial) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0) ;
        vc = update(vc, idxC, idxW, alpha_unr) ;
    end
    if flag_couple
        if correctunCh
            idxC = inputConj(choiceunCh, cnt_trial) ;
            idxW = inputConj(3-choiceunCh, cnt_trial) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0) ;
            vc = update(vc, idxC, idxW, alpha_rew) ;
        else
            idxC = inputConj(3-choiceunCh, cnt_trial) ;
            idxW = inputConj(choiceunCh, cnt_trial) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0) ;
            vc = update(vc, idxC, idxW, alpha_unr) ;
        end
    end
    
    % feature
    if correct
        if cntD==1
            idxC = idx_color(choice) ;
            idxW = idx_color(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_rewColor) ;
        elseif cntD==2
            idxC = idx_shape(choice) ;
            idxW = idx_shape(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_rewShape) ;
        elseif cntD==3
            idxC = idx_pattern(choice) ;
            idxW = idx_pattern(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_rewPattern) ;
        end
    else
        if cntD==1
            idxW = idx_color(choice) ;
            idxC = idx_color(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_unrColor) ;
        elseif cntD==2
            idxW = idx_shape(choice) ;
            idxC = idx_shape(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_unrShape) ;
        elseif cntD==3
            idxW = idx_pattern(choice) ;
            idxC = idx_pattern(3-choice) ;
            [idxW, idxC] = idxcoupleF(idxW, idxC, correct, 0, flag_updatesim) ;
            vf = update(vf, idxC, idxW, alpha_unrPattern) ;
        end
    end
    if flag_couple
        if correctunCh
            if cntD==1
                idxC = idx_color(choiceunCh) ;
                idxW = idx_color(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_rewColor) ;
        elseif cntD==2
                idxC = idx_shape(choiceunCh) ;
                idxW = idx_shape(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_rewShape) ;
        elseif cntD==3
                idxC = idx_pattern(choiceunCh) ;
                idxW = idx_pattern(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_rewPattern) ;
            end
        else
            if cntD==1
                idxW = idx_color(choiceunCh) ;
                idxC = idx_color(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_unrColor) ;
            elseif cntD==2
                idxW = idx_shape(choiceunCh) ;
                idxC = idx_shape(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_unrShape) ;
            elseif cntD==3
                idxW = idx_pattern(choiceunCh) ;
                idxC = idx_pattern(3-choiceunCh) ;
                [idxW, idxC] = idxcoupleF(idxW, idxC, correctunCh, 0, flag_updatesim) ;
                vf = update(vf, idxC, idxW, alpha_unrPattern) ;
            end
        end
    end
    V(1:3,cnt_trial) = vf ;
    V(4:12,cnt_trial) = vc ;
    
end
end


function v = update(v, idxC, idxW, Q)
    if isempty(idxW) && ~isempty(idxC)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
    elseif isempty(idxC) && ~isempty(idxW)
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
    elseif ~isempty(idxW) && ~isempty(idxC)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
    elseif isempty(idxW) && isempty(idxC)
    end
end

function [idxW, idxC] = idxcouple(idxW, idxC, rl2_correct, flag_couple)
    if rl2_correct
        if flag_couple==0
            idxW = [] ;
        end
    else
        if flag_couple==0
            idxC = [] ;
        end
    end
end

function [idxW, idxC] = idxcoupleF(idxW, idxC, rl2_correct, flag_couple, flag_updatesim)
    if rl2_correct
        if flag_couple==0
            idxW = [] ;
        elseif flag_couple==1
            if idxW==idxC                                                  % to avoid potentiating and depressing similar V in coupled cases
                idxW= [] ;
                if ~flag_updatesim
                    idxC = [] ;
                end
            end
        end
    else
        if flag_couple==0
            idxC = [] ;
        elseif flag_couple==1
            if idxW==idxC                                                  % to avoid potentiating and depressing similar V in coupled cases
                idxC= [] ;
                if ~flag_updatesim
                    idxW = [] ;
                end
            end
        end
    end
end