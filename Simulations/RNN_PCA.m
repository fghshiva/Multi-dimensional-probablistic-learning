
clc
clear all
close all

%%

% load('./files/r_vals_PCA.mat')
% 
% idxN         = 1:90 ;   % only on excitatory neurons
% r_valsPCA    = [] ;
% r_valsPCA_av = [] ;
% 
% % concatenate activity/taking the average activity 
% % during the choice period
% for cnt_rep = 1:size(r_vals,3)
%     for cnt_obj = 1:27
%         for t = 1:10
%             r_valsChoice(cnt_obj,t,:) = mean(r_vals((cnt_obj-1+(t-1)*27)*15+[1:15],idxN,cnt_rep),1) ;
%         end
%         r_valsPCA = [r_valsPCA squeeze(r_valsChoice(cnt_obj,:,:))'] ;
%     end
% end
% 
% % concatenate activity/taking the average activity 
% % during the choice period over all repetitions
% for cnt_obj = 1:27
%     for t = 1:10
%         r_valsChoice_av(cnt_obj,t,:) = mean(mean(r_vals((cnt_obj-1+(t-1)*27)*15+[1:15],idxN,:),1),3) ;
%     end
%     r_valsPCA_av = [r_valsPCA_av squeeze(r_valsChoice_av(cnt_obj,:,:))'] ;
% end
% 
% % performing PCA
% [coeff,score,latent,tsquared,explained] = pca(zscore(r_valsPCA,[],2)') ;
% Xcentered = zscore(r_valsPCA_av,[],2)'*coeff ;
% save('./files/r_vals_PCA_choice.mat')

%%

% load('./files/r_vals_PCAtime.mat')
% 
% stps       = [1:2:101] ; % to reduce complexity we subsample
% idxN       = 1:90 ;      % only on excitatory neurons
% r_valsPCA  = [] ;
% r_valsPCA_av = [] ;
% 
% for cnt_rep =  1:size(r_vals,3)
%     for cnt_obj = 1:27
%         r_vals       = r_vals(1:2727,:,:) ;
%         r_valsChoice = r_vals((cnt_obj-1)*101+stps,idxN,cnt_rep) ;
%         r_valsPCA    = [r_valsPCA squeeze(r_valsChoice)'] ;
%     end
% end
% 
% for cnt_obj = 1:27
%     r_valsChoice_av = squeeze(mean(r_vals((cnt_obj-1)*101+stps,idxN,:),3)) ;
%     r_valsPCA_av    = [r_valsPCA_av squeeze(r_valsChoice_av)'] ;
% end
% 
% % performing PCA
% [coeff,score,latent,tsquared,explained] = pca(zscore(r_valsPCA,[],2)') ;
% Xcentered = zscore(r_valsPCA_av,[],2)'*coeff ;
% save('./files/r_vals_PCA_init.mat')

%%

% load('./files/r_vals_PCAtime.mat')
% 
% stps       = [1:2:101] ; % to reduce complexity we subsample
% idxN       = 1:90 ;      % only on excitatory neurons
% r_valsPCA    = [] ;
% r_valsPCA_av = [] ;
% 
% for cnt_rep =  1:size(r_vals,3)
%     for cnt_obj = 1:27
%         r_vals       = r_vals(end-2727:end,:,:) ;
%         r_valsChoice = r_vals((cnt_obj-1)*101+stps,idxN,cnt_rep) ;
%         r_valsPCA    = [r_valsPCA squeeze(r_valsChoice)'] ;
%     end
% end
% 
% for cnt_obj = 1:27
%     r_valsChoice_av = squeeze(mean(r_vals((cnt_obj-1)*101+stps,idxN,:),3)) ;
%     r_valsPCA_av    = [r_valsPCA_av squeeze(r_valsChoice_av)'] ;
% end
% 
% % performing PCA
% [coeff,score,latent,tsquared,explained] = pca(zscore(r_valsPCA,[],2)') ;
% Xcentered = zscore(r_valsPCA_av,[],2)'*coeff ;
% save('./files/r_vals_PCA_end.mat')

%%

close all
prb = [0.92, 0.16, 0.92, 0.75, 0.75, 0.75, 0.43, 0.98, 0.43, ...
       0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , ...
       0.57, 0.02, 0.57, 0.25, 0.25, 0.25, 0.08, 0.84, 0.08] ;
[srt, idxsrt] = sort(unique(prb)) ;

% plotting/coloring based on the reward value
for cnt_obj = 1:27
    idxsrtplot(cnt_obj) = idxsrt(srt==prb(cnt_obj)) ;
end
clrmat = parula(length(idxsrt)) ;

figure(1)
colormap('parula')
hold on
for cnt_obj = 1:27
    xplot = Xcentered((cnt_obj-1)*10+1:cnt_obj*10,[1 2 3]) ;
    if ismember(cnt_obj, [1:18 20 22:25 27])
        plot3(xplot(:,1), xplot(:,2), xplot(:,3), 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    else
        plot3(xplot(:,1), xplot(:,2)*0.4, xplot(:,3), 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    end
    for cnt_plot = 1:1:10
        if ismember(cnt_obj, [1:18 20 22:25 27])
            scatter3(xplot(cnt_plot,1), xplot(cnt_plot,2), xplot(cnt_plot,3), 's', 'SizeData', 8*cnt_plot, ...
                'MarkerFaceColor', clrmat(idxsrtplot(cnt_obj),:), 'MarkerEdgeColor', 'none', 'linewidth', 2)
        else
            scatter3(xplot(cnt_plot,1), xplot(cnt_plot,2)*0.4, xplot(cnt_plot,3), 's', 'SizeData', 8*cnt_plot, ...
                'MarkerFaceColor', clrmat(idxsrtplot(cnt_obj),:), 'MarkerEdgeColor', 'none', 'linewidth', 2)
        end
        
    end
end
view(500,20)
grid on

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'xtick',[-14:7:21],'ytick',[-6:3:9],'ztick',[-6:3:9],'tickdir', 'out','fontsize', 16)
axis([-14, 21, -6, 9, -6, 9])
box off

FigW = 5 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/PCA.pdf')

%%

load('./files/r_vals_PCA_init.mat')
prb = [0.92, 0.16, 0.92, 0.75, 0.75, 0.75, 0.43, 0.98, 0.43, ...
       0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , ...
       0.57, 0.02, 0.57, 0.25, 0.25, 0.25, 0.08, 0.84, 0.08] ;
[srt, idxsrt] = sort(unique(prb)) ;
for cnt_obj = 1:27
    idxsrtplot(cnt_obj) = idxsrt(srt==prb(cnt_obj)) ;
end
clrmat = parula(length(idxsrt)) ;
mrk = {'d', 's', 'v'} ;

figure(1)
colormap('parula')
hold on
for cnt_obj = 1:27
    xplot = Xcentered((cnt_obj-1)*length(stps)+1:cnt_obj*length(stps),[1 2 3]) ;
    cnt_plott(1) = find( abs(51-stps)==min(abs(51-stps)),1 ) ;
    cnt_plott(2) = find( abs(76-stps)==min(abs(76-stps)),1 ) ;
    cnt_plott(3) = find( abs(101-stps)==min(abs(101-stps)),1 ) ;
    idx1 = 1:cnt_plott(1) ;
    idx2 = cnt_plott(1):cnt_plott(2) ;
    idx3 = cnt_plott(2):cnt_plott(3)  ;
    plot3(xplot(idx1,1), xplot(idx1,2), xplot(idx1,3), '--', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    plot3(xplot(idx2,1), xplot(idx2,2), xplot(idx2,3), '-', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    plot3(xplot(idx3,1), xplot(idx3,2), xplot(idx3,3), '--', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    cnt_mrk = 0 ;
    for cnt = [26 51 76]
        cnt_mrk = cnt_mrk+1 ;
        cnt_plot = find( abs(cnt-stps)==min(abs(cnt-stps)),1) ;
        scatter3(xplot(cnt_plot,1), xplot(cnt_plot,2), xplot(cnt_plot,3), mrk{cnt_mrk}, 'SizeData', 50, ...
            'MarkerFaceColor', clrmat(idxsrtplot(cnt_obj),:), 'MarkerEdgeColor', 'none', 'linewidth', 2)
    end
end
view(125,15)
grid on

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'xtick',[-16:8:16],'ytick',[-4:4:8],'ztick',[-6:3:6],'tickdir', 'out','fontsize', 18)
axis([-16, 16, -4, 8, -6, 6])
box off

FigW = 5 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/PCA_init.pdf')

%%

load('./files/r_vals_PCA_end.mat')
prb = [0.92, 0.16, 0.92, 0.75, 0.75, 0.75, 0.43, 0.98, 0.43, ...
       0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , 0.5 , ...
       0.57, 0.02, 0.57, 0.25, 0.25, 0.25, 0.08, 0.84, 0.08] ;
[srt, idxsrt] = sort(unique(prb)) ;
for cnt_obj = 1:27
    idxsrtplot(cnt_obj) = idxsrt(srt==prb(cnt_obj)) ;
end
clrmat = parula(length(idxsrt)) ;
mrk = {'d', 's', 'v'} ;

figure(1)
colormap('parula')
hold on
for cnt_obj = 1:27
    xplot = Xcentered((cnt_obj-1)*length(stps)+1:cnt_obj*length(stps),[1 2 3]) ;
    cnt_plott(1) = find( abs(51-stps)==min(abs(51-stps)),1 ) ;
    cnt_plott(2) = find( abs(76-stps)==min(abs(76-stps)),1 ) ;
    cnt_plott(3) = find( abs(101-stps)==min(abs(101-stps)),1 ) ;
    idx1 = 1:cnt_plott(1) ;
    idx2 = cnt_plott(1):cnt_plott(2) ;
    idx3 = cnt_plott(2):cnt_plott(3)  ;
    plot3(xplot(idx1,1), xplot(idx1,2), xplot(idx1,3), '--', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    plot3(xplot(idx2,1), xplot(idx2,2), xplot(idx2,3), '-', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    plot3(xplot(idx3,1), xplot(idx3,2), xplot(idx3,3), '--', 'color', clrmat(idxsrtplot(cnt_obj),:), 'linewidth', 2)
    cnt_mrk = 0 ;
    for cnt = [26 51 76]
        cnt_mrk = cnt_mrk+1 ;
        cnt_plot = find( abs(cnt-stps)==min(abs(cnt-stps)),1) ;
        scatter3(xplot(cnt_plot,1), xplot(cnt_plot,2), xplot(cnt_plot,3), mrk{cnt_mrk}, 'SizeData', 50, ...
            'MarkerFaceColor', clrmat(idxsrtplot(cnt_obj),:), 'MarkerEdgeColor', 'none', 'linewidth', 2)
    end
end
view(125,15)
grid on

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
set(gca,'xtick',[-20:10:50],'ytick',[-12:6:12],'ztick',[-5:5:10],'tickdir', 'out','fontsize', 18)
axis([-10, 30, -12, 12, -5, 10])
box off

FigW = 5 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/PCA_end.pdf')

%%

figure(3)
colormap('parula')
cb = colorbar() ;
cb.LineWidth = 2 ;
cb.FontSize = 24 ;
set(gca,'xtick',[],'ytick',[],'tickdir', 'out','fontsize', 24, 'linewidth', 2)

FigW = 5 ;
FigH = 4 ;
set(gcf,'units','centimeters')
set(gcf,'position',[10,10,3*FigW,3*FigH],'PaperSize',[FigW FigH],'PaperPosition',...
[0,0,FigW,FigH],'units','centimeters'); 
% print('-dpdf', './figures/clrbr.pdf')
