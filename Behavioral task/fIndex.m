function [idxFtMd idxObjMd idxConjMd idxModel] = fIndex(AIC)
    [~, idxModel]       = min(AIC(:,1:9), [], 2) ;
    idxFtMd             = find(ismember(idxModel, [1 2 3])) ;
    idxObjMd            = find(ismember(idxModel, [4 5 6])) ;
    idxConjMd{1}        = find(ismember(idxModel, [7 8 9])) ;
%     idxConjMd{2}        = find(ismember(idxModel, [10 11 12])) ;
%     idxConjMd{3}        = find(ismember(idxModel, [13 14 15])) ;
        
%     [~, idxModel]       = min(AIC(:, [1 4 7]) , [], 2) ;
%     idxModel            = (idxModel-1)*3 + 1 ;
%     idxFtMd             = find(ismember(idxModel, [1])) ;
%     idxObjMd            = find(ismember(idxModel, [4])) ;
%     idxConjMd{1}        = find(ismember(idxModel, [7])) ;

%     [~, idxModel]       = min(AIC(:, [2 5 8]) , [], 2) ;
%     idxModel            = (idxModel-1)*3 + 2 ;
%     idxFtMd             = find(ismember(idxModel, [2])) ;
%     idxObjMd            = find(ismember(idxModel, [5])) ;
%     idxConjMd{1}        = find(ismember(idxModel, [8])) ;

%     [~, idxModel]       = min(AIC(:, [3 6 9]) , [], 2) ;
%     idxModel            = (idxModel-1)*3 + 3 ;
%     idxFtMd             = find(ismember(idxModel, [3])) ;
%     idxObjMd            = find(ismember(idxModel, [6])) ;
%     idxConjMd{1}        = find(ismember(idxModel, [9])) ;

%     [~, idxModel]       = min(AIC(:,1:9), [], 2) ;
%     idxFtMd             = find(ismember(idxModel, [3])) ;
%     idxObjMd            = find(ismember(idxModel, [6])) ;
%     idxConjMd{1}        = find(ismember(idxModel, [9])) ;
end