function PlotGroupVelocity(GroupToPlot)

% Written by: Loranzie S. Rogers
%
% Will plot fish paths and difference angles for a specified group with
% fish velocity during localizations represented as the path.

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
xPos = -50;
yPos = 60;
radiusRB = 15;
th = 0:pi/50:2*pi;
xCircRB = radiusRB * cos(th) + xPos;
yCircRB = radiusRB * sin(th) + yPos;



% -------------------------------------------------------------------------
%
% Plots the path and difference angle for individual fish. However, the
% velocity of the fish overlays the fish path to the speaker.
%
% -------------------------------------------------------------------------

for i = 1:numel(GroupToPlot)

    xx=[GroupToPlot{1,i}(:,1) GroupToPlot{1,i}(:,1)];
    yy=[GroupToPlot{1,i}(:,2) GroupToPlot{1,i}(:,2)];
    vv = [GroupToPlot{1,i}(:,7) GroupToPlot{1,i}(:,7)];
    zz=zeros(size(xx));


    figure(i); clf;
    t = tiledlayout(2,4);
    setappdata(gcf, 'SubplotDefaultAxesLocation', [0.1, 0.1, 0.85, 1]);

    nexttile([2,2])

    % experimental tank
    plot(xCirc, yCirc,'-k','LineWidth',2)
    axis equal
    yline(0,'--k','LineWidth',2,'Alpha',0.25)
    hold on

    % acoustic monopole
    quiver(x, y, -u, -v,'k')

    % release basket
    plot(xCircRB, yCircRB,'-k','LineWidth',2)

    % speaker
    plot(0,-1,'ok','MarkerFaceColor','k','MarkerSize',8)

    % path with velocity
    hs=surf(xx,yy,zz,vv,'EdgeColor','interp','LineWidth',4);
    colormap('turbo')
    c=colorbar;
    c.Limits = [0 round(max(max(vv)))];
    c.Label.String = 'Velocity (cm/s)';
    c.Label.FontName = 'Arial';
    c.Label.FontSize = 12;
    c.TickDirection = 'out';
    c.LineWidth = 1;
    view(2) %// view(0,90)
    set(gca, 'YDir','reverse')
    ylim([-5 radius+5])
    xlim([-radius-1 radius+1])
    set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',12,'LineWidth',1)
    xlabel('X (cm)','FontName','Arial','FontSize',13)
    ylabel('Y (cm)','FontName','Arial','FontSize',13)
    hold off



    % -------------------------------------------------------------------------
    %
    % Moving mean and sem of angle of error relative to speaker.
    %
    % -------------------------------------------------------------------------

    DistanceToSpeaker=[];
    temp = GroupToPlot{i}(:,5);
    DistanceToSpeaker=[DistanceToSpeaker;temp];

    AngleOfError=[];
    temp = GroupToPlot{i}(:,6);
    AngleOfError=[AngleOfError;temp];

    new = [DistanceToSpeaker,AngleOfError];
    new2 = sortrows(new,1);

    AvgWindow = movmean(new2,5); % use 25 for all except usl, which will use 5
    StdWindow = movstd(new2,5)/sqrt(size(new2,2));

    curve1 = AvgWindow(:,2) + StdWindow(:,2);
    curve2 = AvgWindow(:,2) - StdWindow(:,2);


    % Orientation error plot
    nexttile([2,2])

    hold on
    patch([AvgWindow(:,1); flipud(AvgWindow(:,1))],[curve2; flipud(curve1)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    plot(AvgWindow(:,1),AvgWindow(:,2), 'color', 'r','linewidth',2);
    plot(GroupToPlot{1,i}(:,5),GroupToPlot{1,i}(:,6),'ok','MarkerSize',7,'LineWidth',1)
    set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',12,'XDir','reverse','LineWidth',1)
    ylabel('Orientation error (deg)','FontName','Arial','FontSize',13)
    xlabel('Distance to speaker (cm)','FontName','Arial','FontSize',13)
    ylim([-180 180]);
    yticks([-180:45:180]);

    clear xx yy vv

    t.TileSpacing = 'compact';
    t.Padding = 'compact';

end


end