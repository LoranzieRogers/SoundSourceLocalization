function PlotGroup(GroupToPlot)


% Written by: Loranzie S. Rogers
% 
% Plot of all fish paths and difference angles for a specified experimental 
% group.


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
[x,y] = meshgrid(-120:8:120,0:8:120);
r = sqrt(x.^2 + y.^2); % r in function of (x, y)
r=flipud(r);
theta = atan2(y, x); % theta in function of (x, y)
u = r.*cos(theta); % x component of the vector field
v = r.*sin(theta); % y component of the vector field

% release basket
xPos = -50;
yPos = 60;
radiusRB = 15;
th = 0:pi/50:2*pi;
xCircRB = radiusRB * cos(th) + xPos;
yCircRB = radiusRB * sin(th) + yPos;



% -------------------------------------------------------------------------
% 
% Plots
% 
% -------------------------------------------------------------------------


figure();
t = tiledlayout(3,2);

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
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',10)


nexttile([1,2])
for j = 1:numel(GroupToPlot)
% Difference angle 
plot(GroupToPlot{1,j}(:,5),GroupToPlot{1,j}(:,6),'ok','MarkerSize',7,'LineWidth',1)
hold on
end
ylabel('Difference angle (deg)','FontName','Arial','FontSize',15)
xlabel('Distance to speaker (cm)','FontName','Arial','FontSize',15)
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',10,'XDir','reverse')
ylim([-180 180])
yticks([-180:45:180])

t.TileSpacing = 'compact';
t.Padding = 'compact';

end