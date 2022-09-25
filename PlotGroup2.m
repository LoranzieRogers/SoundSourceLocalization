function PlotGroup2(GroupToPlot)


% Written by: Loranzie S. Rogers
% 
% Plot of all fish paths and difference angles for a specified experimental 
% group.


% -------------------------------------------------------------------------
% 
% Moving mean and sem of angle of error relative to speaker.
% 
% -------------------------------------------------------------------------


DistanceToSpeaker=[];
for i=1:length(GroupToPlot)
   temp = GroupToPlot{i}(:,5);
   DistanceToSpeaker=[DistanceToSpeaker;temp];
end

AngleOfError=[];
for i=1:length(GroupToPlot)
   temp = GroupToPlot{i}(:,6);
   AngleOfError=[AngleOfError;temp];
end


new = [DistanceToSpeaker,AngleOfError];
new2 = sortrows(new,1);

AvgWindow = movmean(new2,5); % use 25 for all except usl, which will use 5
StdWindow = movstd(new2,25)/sqrt(size(new2,2));
% StdWindow = movstd(new2,30);

curve1 = AvgWindow(:,2) + StdWindow(:,2);
curve2 = AvgWindow(:,2) - StdWindow(:,2);


% -------------------------------------------------------------------------
% 
% Information about the tank and monopole speaker.
% 
% -------------------------------------------------------------------------

% experimental tank
radius=176;
ArenaTheta = linspace(2*pi/2, 0, 100);  % <-- left half of circle
xCirc = radius * cos(ArenaTheta);
yCirc = radius * sin(ArenaTheta);

% acoustic monopole
[x,y] = meshgrid(-200:15:200,0:15:280);
r = sqrt(x.^2 + y.^2); % r in function of (x, y)
r=flipud(r);
theta = atan2(y, x); % theta in function of (x, y)
u = r.*cos(theta); % x component of the vector field
v = r.*sin(theta); % y component of the vector field

% release basket
xPos = -44;
yPos = 64;
radiusRB = 12;
th = 0:pi/50:2*pi;
xCircRB = radiusRB * cos(th) + xPos;
yCircRB = radiusRB * sin(th) + yPos;



% -------------------------------------------------------------------------
% 
% Plots
% 
% -------------------------------------------------------------------------


figure();
t = tiledlayout(2,4);

nexttile([2,2])
% experimental tank
plot(xCirc, yCirc,'-k','LineWidth',2)
axis equal
yline(0,'--k','LineWidth',2,'Alpha',0.25)
hold on

% acoustic monopole
quiver(x, y, -u, -v,'k')

% speaker
plot(0,-1,'ok','MarkerFaceColor','k','MarkerSize',8)

% release basket
plot(xCircRB, yCircRB,'-k',"LineWidth",2)


for i = 1:numel(GroupToPlot)
% path of fish to speaker
plot(GroupToPlot{1,i}(:,1),GroupToPlot{1,i}(:,2),'-','LineWidth',1.5)
hold on
end
set(gca, 'YDir','reverse')
ylim([-5 radius+5])
xlim([-radius radius])
xlabel('X (cm)','FontName','Arial','FontSize',15)
ylabel('Y (cm)','FontName','Arial','FontSize',15)
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',13.5,'LineWidth',1)
yticks(0:20:180);


nexttile([2,2])

hold on;
patch([AvgWindow(:,1); flipud(AvgWindow(:,1))],[curve2; flipud(curve1)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
plot(AvgWindow(:,1),AvgWindow(:,2), 'color', 'r','linewidth',2);
plot(DistanceToSpeaker,AngleOfError,'ok','MarkerSize',7,'LineWidth',1)
hold on
ylabel('Difference angle (deg)','FontName','Arial','FontSize',15)
xlabel('Distance to speaker (cm)','FontName','Arial','FontSize',15)
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',13.5,'XDir','reverse','LineWidth',1)
xlim([0 100]);
xticks(0:10:100);
ylim([-180 180])
yticks(-180:45:180)

t.TileSpacing = 'compact';
t.Padding = 'compact';

end