function loglikehood = fMLchoicefit_RL2objdecay(xpar, sesdata)
%
% DESCRIPTION: fits data to RL(2)obj model using ML method
%
% INPUT: 
% sesdata structure which includes input, experiment and behavioral data
%
% OUTPUT:
% fitted parametres

loglikehood = 0 ;
NparamBasic = 3 ;

BiasL   = xpar(1) ;
mag     = xpar(2) ;

xpar([NparamBasic:NparamBasic+sesdata.Nalpha])=1./(1+exp(-(xpar([NparamBasic:NparamBasic+sesdata.Nalpha]))./sesdata.sig) ) ;
decay = xpar(3) ;
alpha_rew = xpar([NparamBasic+1]) ;
if sesdata.flagUnr==1
    alpha_unr = xpar([NparamBasic+2]) ;
else
    alpha_unr = alpha_rew ;
end

inputTarget         = sesdata.input.inputTarget ;
correcttrials       = sesdata.results.reward ;
choicetrials        = sesdata.results.choice ;
flag_couple         = sesdata.flag_couple ;
ntrials             = length(choicetrials) ;
inputRewards        = sesdata.input.inputReward ;

v = (0.5*ones(27,1)) ; 
for cnt_trial=1:ntrials
    
    correct = correcttrials(cnt_trial) ;
    choice = choicetrials(cnt_trial) ; 
    correctunCh = inputRewards(3-choice, cnt_trial) ;
    choiceunCh = 3-choice ;
    
    pChoiceR = 1./(1+exp(-( mag*(v(inputTarget(2, cnt_trial))-v(inputTarget(1, cnt_trial))) + BiasL ) )) ;
    pChoiceL = 1-pChoiceR ;
    if cnt_trial >= 1  
        if choice == 2 
            loglikehood = loglikehood - log(pChoiceR) ;
        else 
            loglikehood = loglikehood - log(pChoiceL) ; 
        end                      
    end
    
    if correct
        idxC = inputTarget(choice, cnt_trial) ;
        idxW = inputTarget(3-choice, cnt_trial) ;
        v = decayV(v, find([1:27]~=inputTarget(choice, cnt_trial)), decay) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0) ;
        v = update(v, idxC, idxW, alpha_rew) ;
    else
        idxC = inputTarget(3-choice, cnt_trial) ;
        idxW = inputTarget(choice, cnt_trial) ;
        v = decayV(v, find([1:27]~=inputTarget(choice, cnt_trial)), decay) ;
        [idxW, idxC] = idxcouple(idxW, idxC, correct, 0) ;
        v = update(v, idxC, idxW, alpha_unr) ;
    end
    if flag_couple
        if correctunCh
            idxC = inputTarget(choiceunCh, cnt_trial) ;
            idxW = inputTarget(3-choiceunCh, cnt_trial) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0) ;
            v = update(v, idxC, idxW, alpha_rew) ;
        else
            idxC = inputTarget(3-choiceunCh, cnt_trial) ;
            idxW = inputTarget(choiceunCh, cnt_trial) ;
            [idxW, idxC] = idxcouple(idxW, idxC, correctunCh, 0) ;
            v = update(v, idxC, idxW, alpha_unr) ;
        end
    end
    V(:,cnt_trial) = v ;
    
end
end

function v = decayV(v, unCh, decay)
	v(unCh) = v(unCh) - (v(unCh)-0.5)*(decay) ;
end

function v = update(v, idxC, idxW, Q)
    if isempty(idxW)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
    elseif isempty(idxC)
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
    elseif ~isempty(idxW) && ~isempty(idxC)
        v(idxC) = v(idxC) + (1-v(idxC)).*Q ;
        v(idxW) = v(idxW) - (v(idxW).*Q) ;
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