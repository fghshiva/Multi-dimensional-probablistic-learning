function loglikehood = fMLchoicefit_RL2ft(xpar, sesdata)
%
% DESCRIPTION: fits data to RL(2) model using ML method
%
% INPUT: 
% sesdata structure which includes input, experiment and behavioral data
%
% OUTPUT:
% fitted parametres

loglikehood = 0 ;
NparamBasic = 4 ;

xpar(2:4)  = abs(xpar(2:4)) ;

BiasL      = xpar(1) ;
magColor   = xpar(2) ;
magShape   = xpar(3) ;
magPattern = xpar(4) ;

xpar([NparamBasic+1:NparamBasic+sesdata.Nalpha])=1./(1+exp(-(xpar([NparamBasic+1:NparamBasic+sesdata.Nalpha]))./sesdata.sig) ) ;
alpha_rewColor      = xpar([NparamBasic+1]) ;
alpha_rewShape      = xpar([NparamBasic+1]) ;
alpha_rewPattern    = xpar([NparamBasic+1]) ; 

alpha_unrColor      = xpar([NparamBasic+2]) ;
alpha_unrShape      = xpar([NparamBasic+2]) ;
alpha_unrPattern    = xpar([NparamBasic+2]) ;

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

v = (0.5*ones(9,1)) ; 
for cnt_trial=1:ntrials
    
    correct = correcttrials(cnt_trial) ;
    choice  = choicetrials(cnt_trial) ; 
    correctunCh = inputRewards(3-choice, cnt_trial) ;
    choiceunCh  = 3-choice ;
    
    idx_shape(2)   = shapeMap(inputTarget(2, cnt_trial)) ;
    idx_color(2)   = colorMap(inputTarget(2, cnt_trial)) ;
    idx_pattern(2) = patternMap(inputTarget(2, cnt_trial)) ;
    idx_shape(1)   = shapeMap(inputTarget(1, cnt_trial)) ;
    idx_color(1)   = colorMap(inputTarget(1, cnt_trial)) ;
    idx_pattern(1) = patternMap(inputTarget(1, cnt_trial)) ;
    
    pChoiceR = 1./(1+exp(-( magShape*(v(idx_shape(2))-v(idx_shape(1))) + magColor*(v(idx_color(2))-v(idx_color(1))) + ...
        magPattern*(v(idx_pattern(2))-v(idx_pattern(1))) + BiasL ) )) ;
    pChoiceL = 1-pChoiceR ;
    if cnt_trial >= 1  
        if choice == 2 
            loglikehood = loglikehood - log(pChoiceR) ;
        else
            loglikehood = loglikehood - log(pChoiceL) ; 
        end                      
    end
    
    if correct
        idxC = idx_color(choice) ;
        idxW = idx_color(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_rewColor) ;

        idxC = idx_shape(choice) ;
        idxW = idx_shape(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_rewShape) ;
        
        idxC = idx_pattern(choice) ;
        idxW = idx_pattern(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_rewPattern) ;
    else
        idxW = idx_color(choice) ;
        idxC = idx_color(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_unrColor) ;

        idxW = idx_shape(choice) ;
        idxC = idx_shape(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_unrShape) ;
        
        idxW = idx_pattern(choice) ;
        idxC = idx_pattern(3-choice) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0, flag_updatesim) ;
        v = update(v, idxC, idxW, alpha_unrPattern) ;
    end
    if flag_couple
        if correctunCh
            idxC = idx_color(choiceunCh) ;
            idxW = idx_color(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_rewColor) ;

            idxC = idx_shape(choiceunCh) ;
            idxW = idx_shape(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_rewShape) ;
            
            idxC = idx_pattern(choiceunCh) ;
            idxW = idx_pattern(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_rewPattern) ;
        else
            idxW = idx_color(choiceunCh) ;
            idxC = idx_color(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_unrColor) ;

            idxW = idx_shape(choiceunCh) ;
            idxC = idx_shape(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_unrShape) ;
            
            idxW = idx_pattern(choiceunCh) ;
            idxC = idx_pattern(3-choiceunCh) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0, flag_updatesim) ;
            v = update(v, idxC, idxW, alpha_unrPattern) ;
        end
    end
    V(:,cnt_trial) = v ;
    
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

function [idxW, idxC] = idxcouple(idxW, idxC, rl2_correct, flag_couple, flag_updatesim)
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