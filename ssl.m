%% Plainfin midshipman sound source localization behavioral analysis
% 
% 
% Written by: Loranzie S. Rogers


%% Clean up environment and set directory

clear; close all; clc;                                                  % clear and close everything
cd('/Users/loranzierogers/Documents/UW/projects/SoundSourceLocalization');

%% Read in all files and extract file names

% location of localization tracks
FileLocation = '/Users/loranzierogers/Documents/UW/projects/SoundSourceLocalization';

% load file info
files = dir(fullfile(FileLocation,'PN*.csv'));

% extract and store file names for quick access
FishID = {files.name};

%% Calibrate data

% -------------------------------------------------------------------------
% 
% Here, we initialize the structure 'track' and fill with all file data.
% 
% tracks: stucture that contains individual 2-D track information.
% 
% Each field has five columns, which contain ID, point #, x pos (px),
% y pos (px), and frame # analyzed in mTrackj. 
% 
% ID: speaker (1) and fish (2)
% 
% For each field row 1 references information related to the speaker, while
% all other rows represent track of fish.
% 
% -------------------------------------------------------------------------

tracks=cell(1,length(files));
for i = 1:numel(files)
    tracks{i} = readmatrix(files(i).name);
end

% -------------------------------------------------------------------------
% 
% Set speaker as the origin for all positions such that we have track
% position information (x,y) relative to the speaker, which is the end
% point for all successful trials.
% 
% CalTracks: stuct with individual track information after setting
% speaker as global origin.
% 
% Each field has two columns, which contain x pos (px) and y pos (px). 
% Once speaker (s) has been set as the origin, you will have the following
% coordinate system:
% 
%           |
%   (-x,-y) | (x,-y)
%   --------s--------
%   (-x,y)  | (x,y)
%           |
% 
% -------------------------------------------------------------------------

CalTracks = cell(1,numel(files));
for i = 1:numel(files)
%         a = abs(tracks{1,j}(2:end,3:4)- tracks{1,j}(1,3:4));
    CorrectedPosition = tracks{1,i}(2:end,3:4) - tracks{1,i}(1,3:4);        
    CalTracks{1,i} = CorrectedPosition;

end

% add corrected frame # to c3 of each field relative to 1st frame analyzed
for i = 1:numel(files)
    CalTracks{1,i}(:,3) = tracks{1,i}(2:end,5) - tracks{1,i}(2,5);
end

% frame rate 30 Hz
fr = 1/30;

% add time to c4
for i = 1:numel(files)
    CalTracks{1,i}(:,4) = CalTracks{1,i}(:,3) .* fr;
end

% collect time to reach speaker for each trial
TimeToSpeaker = zeros(size(CalTracks));
for i = 1:numel(CalTracks)
    TimeToSpeaker(1,i) = CalTracks{1,i}(end,4);
end

% average time to speaker
AvgT = mean(TimeToSpeaker,'omitnan');


% coversion from pixels to cm
cal = readmatrix('cal.csv')';

% Convert x pos (px) and y pos (px) to cm use cal factor
for i = 1:numel(files)
    CalTracks{1,i}(:,1) = CalTracks{1,i}(:,1) ./ cal(1,i);
    CalTracks{1,i}(:,2) = CalTracks{1,i}(:,2) ./ cal(1,i);
end


% distance (cm) traveled along path to speaker
DistanceTraveled = zeros(size(CalTracks));
for i = 1:numel(files)
    DistanceTraveled(1,i) = sum(hypot(diff(CalTracks{1,i}(:,1)), diff(CalTracks{1,i}(:,2))));
end



%% Perform analyses

% -------------------------------------------------------------------------
% 
% Determine fish distance to speaker. 
% 
% Distnace between two points is determined using pythagorean theorem.
% 
% The distance (d) between two points is as follows:
% 
% d = sqrt((x2-x1).^2 + (y2-y1).^2)
% 
% Since the position of the speaker (x2,y2) is (0,0) we are left with only
% the position of the fish (x1,y1) and are equation is now: 
% 
% d = sqrt(x1.^2 = y2.^2)
% 
% -------------------------------------------------------------------------

for i = 1:numel(files)
    for j = 1:numel(CalTracks{1,i}(:,1))
        CalTracks{1,i}(j,5) = sqrt(CalTracks{1,i}(j,1).^2 + CalTracks{1,i}(j,2).^2);
    end
end


% -------------------------------------------------------------------------
% 
% Angle relative to sound source.
% 
% Here, we determine the difference angles of the fish bearing relative to 
% the sound source similar to Zeddies et al., 2010, 2012. 
% 
% The bearing angle (theta) is defined as the line segment from the
% position of the speaker (sx,sy) to each point along the fish path 
% (px,py). 
% 
% 
% -------------------------------------------------------------------------

% for i = 1:numel(files)
%     for j = 1:numel(CalTracks{1,i}(:,1))
%         CalTracks{1,i}(j,6) = -atan2(0 - CalTracks{1,i}(j,2), 0 - CalTracks{1,i}(j,1)) * 180 / pi;  % needs to be negative so it is relative to release site. wrapeto180()
%     end
% end


for i = 1:numel(files)
    dxy = diff([CalTracks{1,i}(:,1),CalTracks{1,i}(:,2)]);
    dxy = dxy .* diff(CalTracks{1,i}(:,4));
    theta1 = atan2(CalTracks{1,i}(:,2), CalTracks{1,i}(:,1)) * 180/pi;
    theta2 = atan2(dxy(:,2),dxy(:,1)) * 180/pi;
    CalTracks{1,i}(:,6) = wrapTo180([0; theta1(2:end,:)-(theta2+180)]);

    clear dxy theta1 theta2
end

% -------------------------------------------------------------------------
% 
% Determine the magnitude of fish instantaneous velocity |v|. Here,
% 
% |v| = dr/dt
% |v| = sqrt((delta x^2) + (delta y^2)) / delta t
% |v| = sqrt(((xj - xi)^2) + ((yj - yi)^2)) / (tj - ti)
% 
% where, i equals the position of the fish at a given time and j = i+1 
% 
% Let's also determine average velocity for each of these trials.
% 
% -------------------------------------------------------------------------

% instantaneous velocity
for i = 1:numel(files)

    dvx = zeros(length(CalTracks{1,i}(:,4))-1,1);
    dvy = zeros(length(CalTracks{1,i}(:,4))-1,1);
    dt = zeros(length(CalTracks{1,i}(:,4))-1,1);
    v = zeros(length(CalTracks{1,i}(:,4)),1);
    d = zeros(length(CalTracks{1,i}(:,4))-1,1);

    for j = 1:numel(CalTracks{1,i}(:,4))-1
        dvx(j,1) = CalTracks{1,i}(j+1,1) - CalTracks{1,i}(j,1);
        dvy(j,1) = CalTracks{1,i}(j+1,2) - CalTracks{1,i}(j,2);
        dt(j,1) = CalTracks{1,i}(j+1,4) - CalTracks{1,i}(j,4);
        v(j+1,1) = sqrt( dvx(j,:).^2 + dvy(j,:).^2 ) ./ dt(j,:);
        d(j,1) = sum(hypot(diff(CalTracks{1,i}(j,1)), diff(CalTracks{1,i}(j,2))));
    end

    CalTracks{1,i}(:,7) = v;
    clear dvx dvy dt v

end

% average velocity
AvgV = zeros(size(CalTracks));
for i = 1:numel(CalTracks)
    temp = CalTracks{1,i}(:,7);
    temp(temp == 0) = NaN;
    AvgV(1,i) = mean(temp,'omitnan');
end

clear temp

% -------------------------------------------------------------------------
% 
% Now, let's add instantaneous acceleration |a|. Here, 
% 
% |a| = dv/dt
% |a| = delta v / delta t
% |a| = (vj - vi) / (tj - ti)
% 
% where, i equals the position of the fish at a given time and j = i+1 
% 
% Again, let's also determine average acceleration for all trials.
% 
% -------------------------------------------------------------------------

% instantaneous acceleration
for i = 1:numel(files)

    dvx = zeros(length(CalTracks{1,i}(:,4))-1,1);
    dvy = zeros(length(CalTracks{1,i}(:,4))-1,1);
    dt = zeros(length(CalTracks{1,i}(:,4))-1,1);
    v = zeros(length(CalTracks{1,i}(:,4)),1);
    d = zeros(length(CalTracks{1,i}(:,4))-1,1);

    for j = 1:numel(CalTracks{1,i}(:,4))-1
        dvx(j,1) = CalTracks{1,i}(j+1,1) - CalTracks{1,i}(j,1);
        dvy(j,1) = CalTracks{1,i}(j+1,2) - CalTracks{1,i}(j,2);
        dt(j,1) = CalTracks{1,i}(j+1,4) - CalTracks{1,i}(j,4);
        a(j+1,1) = sqrt( dvx(j,:).^2 + dvy(j,:).^2 ) ./ (dt(j,:).^2);
    end

    CalTracks{1,i}(:,8) = a;
    clear dvx dvy dt a

end


% average acceleration
AvgA = zeros(size(CalTracks));
for i = 1:numel(CalTracks)
    temp = CalTracks{1,i}(:,8);
    temp(temp == 0) = NaN;
    AvgA(1,i) = mean(temp,'omitnan');
end



%% Partition tracks into experimental groups

% -------------------------------------------------------------------------
% 
% Using FishID, we are going to partition each cell from CalTracks into
% their corresponding experimental group. This will allow us to plot
% individual groups together more easily. Experimental groups are as
% follows:
% 
% 
% control = Control (n = 9,31) (# of postive responses)
% bs = Bilateral saccule removal (n = 3)
% bul = Bilateral utricle and lagena removal (n = 2)
% ul = Unilateral lagena removal (n = 7)
% bl = Bilateral lagena removal (n = 3)
% usl = Unilateral saccule and lagena removal (n = 1)
% usul = Unilateral saccule, utricle, and lagena removal (n = 2)
% 
% -------------------------------------------------------------------------

control = CalTracks(1,[1,2,3,5,6,9]);
bs = CalTracks(1,[4,7,8]);
ul = CalTracks(1,[10,12,14,16]);
bl = CalTracks(1,[11,13,15]);
usl = CalTracks(1,17);
usul = CalTracks(1,[18,19]);

% CONTINUE TO ADD ONCE ALL DATA HAS BEEN COLLEECTED....


%% Plots

% -------------------------------------------------------------------------
% 
% There are a few options that we can choose from in order to obtain our
% desired plot outputs for individual experimental groups - each of which 
% have been written as a function to make the code less cluttered. 
% Plotting functions are as follows:
% 
% PlotGroup: Plot of all fish paths and difference angles for a specified 
% experimental group.
% 
% PlotGroupIndividuals: Plots individual fish paths and difference angles 
% for a specified group.
% 
% PlotGroupVelocity: Will plot fish paths and difference angles for a 
% specified group with fish velocity during localizations represented as 
% the path.
% 
% -------------------------------------------------------------------------

PlotGroup(control)
PlotGroup2(control)
PlotGroupIndividuals(control)




%% Further analysis of fish paths to speaker

% concatenate data
DTVA = [DistanceTraveled; TimeToSpeaker; AvgV; AvgA];

% Plot path stats - see function for additional details
PlotPathStats(DTVA)

%% Plot proportion of positive and negative responses

% control, ul, bs, bl, usul, bul, usl, us, bsl
ProportionOfPosResponses = [0.2727 0.2593 0.1250 0.12 0.0833 0.0645 0.0417 0 0];
ProportionOfNegResponse = 1-ProportionOfPosResponses;
ProportionOfResponses = [ProportionOfPosResponses;ProportionOfNegResponse];
x_ax = 1:9;

figure()
b = barh(x_ax, ProportionOfResponses,'stacked','FaceColor','flat');
b(1).CData = [0 0 0];
b(2).CData = [1 1 1]*0.8;
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',13.5,'LineWidth',1,'YDir','reverse')
xlabel('Proportion of Responses','FontName','Arial','FontSize',15)
ylim([0.25 9.75])
yticklabels({'Control','Unilateral lagena','Bilateral saccule','Bilateral lagena','Unilateral ear','Bilateral utricle and lagena','Unilateral saccule and lagena','Unilateral sacuule','Bilateral saccule and lagena'});
xticks(0:0.2:1)


%% EXTRA CODE - Tracking position through time

% This code is a play off of a DeepLabCut figure

plot3(ul{1,2}(:,4),ul{1,2}(:,1),ul{1,2}(:,2))
grid on
pbaspect([4 1 1])
% view(3)
set(gca,'Ydir','reverse','Zdir','reverse')
ylim([-180 180])
yticks(-180:45:180)
ylabel('X (cm)')
zlim([0 180])
zticks(0:20:180)
zlabel('Y (cm)')
xlabel('Time (s)')

%%

% imagesc

for i = 1:numel(control)
% path of fish to speaker
plot(control{1,i}(:,5),control{1,i}(:,7),'-','LineWidth',1.5)
hold on
end


for i = 1:numel(control)
% path of fish to speaker
imagesc(control{1,i}(:,7))
hold on
end

%%
velocity = readmatrix("velocity2.csv");
velocity = velocity';
vel2 = velocity(:,181:end);

smtlb = sgolayfilt(velocity,1,3);
NormVelocity = velocity ./ max(velocity,[],2);

Norm = smtlb ./ max(smtlb,[],2);

figure(1); clf;
imagesc(Norm)
colormap(fire)
c= colorbar;
% xlim([170 192])

%%

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
radiusRB = 12;
th = 0:pi/50:2*pi;
xCircRB = radiusRB * cos(th) + xPos;
yCircRB = radiusRB * sin(th) + yPos;

plot(xCirc, yCirc,'-k','LineWidth',2)
axis equal
yline(0,'--k','LineWidth',2,'Alpha',0.25)
hold on

% acoustic monopole
quiver(x, y, -u, -v,'k')

% speaker
plot(0,-1,'ok','MarkerFaceColor','k','MarkerSize',8)

% release basket
plot(xCircRB, yCircRB,'-k','LineWidth',2)

% path of fish to speaker
a = plot(control{1,1}(:,1),control{1,1}(:,2),'-','LineWidth',2)
hold on
set(gca, 'YDir','reverse')
ylim([-5 radius+5])
xlim([-radius radius])
set(gca,'Box','off','TickDir','out','FontName','Arial','FontSize',12,'LineWidth',1)
xlabel('X (cm)','FontName','Arial','FontSize',13)
ylabel('Y (cm)','FontName','Arial','FontSize',13)
set(a.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

%%

newDistance = 0:1:80;
OriginalVelocity

NewVelocity = interp1(OriginalDistance, oldYvals, newDistance);
