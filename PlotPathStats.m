function PlotPathStats(DTVA)

% -------------------------------------------------------------------------
% 
% PlotPathStats(): Function that takes in the array DTVA, which represent
% the distance traveled (cm), time to speaker (s), velocity (cm/s), and 
% acceleration (cm/s^2) of fish that exhibited positive responses (i.e., 
% localized underwater sound source) from each treatment group. Code will
% subset DTVA into treatment groups based on manually inputs, which agree
% with the appropriate FishID, and then will compute the mean +/- sem of
% each group for all variables of interest. Short hand labels confer with
% those used in main code. 
% 
% Written by: Loranzie S. Rogers
% 
% -------------------------------------------------------------------------


% subset data based on treatment - these are manually done based on input
controlDTVA = DTVA(:,[1,2,3,5,6,9]);
bsDTVA = DTVA(:,[4,7,8]);
ulDTVA = DTVA(:,[10,12,14,16]);
blDTVA = DTVA(:,[11,13,15]);
uslDTVA = DTVA(:,17);
usulDTVA = DTVA(:,[18,19]);

% mean +/- sem
controlMean = mean(controlDTVA,2); controlSEM = std(controlDTVA,[],2)/sqrt(size(controlDTVA,2));
bsMean = mean(bsDTVA,2); bsSEM = std(bsDTVA,[],2)/sqrt(size(bsDTVA,2));
ulMean = mean(ulDTVA,2); ulSEM = std(ulDTVA,[],2)/sqrt(size(ulDTVA,2));
blMean = mean(blDTVA,2); blSEM = std(blDTVA,[],2)/sqrt(size(blDTVA,2));
uslMean = mean(uslDTVA,2); uslSEM = std(uslDTVA,[],2)/sqrt(size(uslDTVA,2));
usulMean = mean(usulDTVA,2); usulSEM = std(usulDTVA,[],2)/sqrt(size(usulDTVA,2));
TreatmentMean = [controlMean bsMean ulMean blMean uslMean usulMean];
TreatmentSEM = [controlSEM bsSEM ulSEM blSEM uslSEM usulSEM];

% create x-axis for plotting
Xcontrol = ones(size(controlDTVA,2),1);
Xbs = 2.*ones(size(bsDTVA,2),1);
Xul = 3.*ones(size(ulDTVA,2),1);
Xbl = 4.*ones(size(blDTVA,2),1);
Xusl = 5.*ones(size(uslDTVA,2),1);
Xusul = 6.*ones(size(usulDTVA,2),1);
x_ax = 1:6; 

% -------------------------------------------------------------------------
% 
% Plot - 2x2 figure that plots distance traveled (cm), time to speaker (s),
% velocity (cm/s), and acceleration (cm/s^2) of fish that exhibited
% positive responses (i.e., localized underwater sound source) from each
% treatment group. Will plot individual data points as black dots and the
% mean +/- sem for each group in red.
% 
% -------------------------------------------------------------------------

figure(); clf;
t = tiledlayout(2,2);

% Distance (cm) traveled during underwater sound source localization 
nexttile()
hold on
scatter(Xcontrol,controlDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbs,bsDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xul,ulDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbl,blDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusl,uslDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusul,usulDTVA(1,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
errorbar(x_ax,TreatmentMean(1,:),TreatmentSEM(1,:),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8,'LineWidth',1.25,'CapSize',0);
set(gca,'box','off','tickdir','out','TickLength',[0.02 0.02],'FontName','Arial','fontsize',13,'LineWidth',1.25);
xlim([0.5,6.5]);
ylim([0 400]);
yticks(0:50:400);
xticks([1:1:6]);
% xticklabels({'Control','Bilateral saccule','Unilateral lagena','Bilateral lagena','Unilateral sul'});
xticklabels({});
ylabel ('Distance (cm)','FontName','Arial','fontsize',15);
hold off

% time to localize speaker (s)
nexttile()
hold on
scatter(Xcontrol,controlDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbs,bsDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xul,ulDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbl,blDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusl,uslDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusul,usulDTVA(2,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
errorbar(x_ax,TreatmentMean(2,:),TreatmentSEM(2,:),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8,'LineWidth',1.25,'CapSize',0);
set(gca,'box','off','tickdir','out','TickLength',[0.02 0.02],'FontName','Arial','fontsize',13,'LineWidth',1.25);
xlim([0.5,6.5]);
ylim([0 40]);
yticks(0:5:40);
xticks([1:1:6]);
% xticklabels({'Control','Bilateral saccule','Unilateral lagena','Bilateral lagena','Unilateral sul'});
xticklabels({});
ylabel ('Time (s)','FontName','Arial','fontsize',15);
hold off

% Velocity of fish (cm/s)
nexttile()
hold on
scatter(Xcontrol,controlDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbs,bsDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xul,ulDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbl,blDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusl,uslDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusul,usulDTVA(3,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
errorbar(x_ax,TreatmentMean(3,:),TreatmentSEM(3,:),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8,'LineWidth',1.25,'CapSize',0);
set(gca,'box','off','tickdir','out','TickLength',[0.02 0.02],'FontName','Arial','fontsize',13,'LineWidth',1.25);
xlim([0.5,6.5]);
ylim([0 30]);
yticks(0:5:30);
xticks([1:1:6]);
xticklabels({'Control','Bilateral saccule','Unilateral lagena','Bilateral lagena','Unilateral saccule and lagena','Unilateral ear'});
ylabel ('Velocity (cm/s)','FontName','Arial','fontsize',15);
hold off

% acceleration of fish (cm/s^2)
nexttile()
hold on
scatter(Xcontrol,controlDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbs,bsDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xul,ulDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xbl,blDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusl,uslDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
scatter(Xusul,usulDTVA(4,:)',150,'.k','MarkerFaceColor','k','jitter','on','jitterAmount',0.3);
errorbar(x_ax,TreatmentMean(4,:),TreatmentSEM(4,:),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8,'LineWidth',1.25,'CapSize',0);
set(gca,'box','off','tickdir','out','TickLength',[0.02 0.02],'FontName','Arial','fontsize',13,'LineWidth',1.25);
xlim([0.5,6.5]);
ylim([0 200]);
yticks(0:50:200);
xticks([1:1:6]);
xticklabels({'Control','Bilateral saccule','Unilateral lagena','Bilateral lagena','Unilateral saccule and lagena','Unilateral ear'});
ylabel ('Acceleration (cm/s^2)','FontName','Arial','fontsize',15);
hold off

t.TileSpacing = 'compact';

end