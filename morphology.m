%% morphology.m
% This scripts takes the .csv file from ImageJ morphology analysis and
% separates values per fibril.

% For this script to work, make sure to measure a very small distance prior
% to every new fibril (treshold < 1.5 nm)

% Written by Elliot J. Crooks

%% A. Inputs: Read file & Threshold
clear; close; clc;

pth(:)='~/Desktop/A2504_6/'; % /path/to/folder/
job(:)='Results.csv'; % .csv file
T= readtable(strcat(pth,job));  % open CSV and store as table
threshold=1.5;      

%% B. Segregate by fibril
columnLength=find(strcmp(T.Properties.VariableNames,'Length'));   % Find which column in the table lengths values are stored

% Find fibril start point

totalFibrils=1;   
for a=1:height(T)   % Iterating through all rows of table
    if T{a,columnLength}<= threshold   % Checking whether this line notes a new fibril starting
        helixStart(totalFibrils)=a; totalFibrils=totalFibrils+1; end    % Storing which lines each new fibril starts at
end

W=T;    % Making a working table

% Removing rows containing values below threshold = 1.5
for a=totalFibrils-1:-1:1; W(helixStart(a),:)=[]; end

% Taking into account row removal in tracking helix start
for a=1:totalFibrils-1; helixStart(a)=helixStart(a)-(a-1); end

% If user forgot to note start of first fibril, this line makes a note of
% it and removes 1 from total fibril count
if helixStart(1)~=1; helixStart=[1, helixStart]; else totalFibrils=totalFibrils-1; end

% If user missed one measurement (either a pitch or width measurement), an
% error message appears:
% if rem(height(W),2)~=0
%     error('One of the measurements is missing its paired value. Go into table T and delete the unpaired row.'); end


% Making a table to store all values without segregating by fibrils
for a=1:2:height(W)
    I(a,1)=W{a,columnLength}; 
    I(a,2)=W{a+1,columnLength}; 
end

for a=height(I)-1:-2:2; I(a,:)=[]; end  % Removing empty rows
I=table(I);I=splitvars(I,'I');
% Tracking helix start having separated width and length values
helixStartNew=(helixStart+1)/2;

% Separating by fibril
for a=1:totalFibrils
    if a~=totalFibrils
        for b=helixStartNew(a):(helixStartNew(a+1)-1)
            F(b-(helixStartNew(a)-1),2*a-1)=I{b, 1};
            F(b-(helixStartNew(a)-1),2*a)=I{b, 2}; end; end
    if a==totalFibrils
         for b=helixStartNew(a):height(I)
            F(b-(helixStartNew(a)-1),2*a-1)=I{b, 1};
            F(b-(helixStartNew(a)-1),2*a)=I{b, 2}; end; end  
end
F=table(F); F=splitvars(F,'F');    
    
%% Calculations

% Finding which row of F the average & stdev will be placed
rowAverage=height(F)+1;
rowStdev=height(F)+2;

% Average and standard deviation of fibril pitch
for a=1:width(F)
    F{rowAverage,a}=mean(nonzeros(F{1:(rowAverage-1),a}));   % average per fibril
    if length(nonzeros(F{1:(rowAverage-1),a}))==1
        F{rowStdev,a}=0; 
    else F{rowStdev,a}=std(nonzeros(F{1:(rowAverage-1),a})); end    % standard deviation per fibril
end

% Average pitch and stddev per fibril: columns 1 and 2 in Fs
counter=1;
for a=1:2:width(F)
    Fs(counter,1)=F{rowAverage,a}; Fs(counter,2)=F{rowStdev,a}; counter=counter+1; end

% Average width and stddev per fibril: columns 3 and 4 in Fs
for a=2:2:width(F)
    Fs(a/2,3)=F{rowAverage,a}; Fs(a/2,4)=F{rowStdev,a}; end

%% Plots

% Initial plots without fibril segregation

f1=figure;
s1=scatter(I{:,1},I{:,2},'filled');
title(['Morphology per measurement. N=', num2str(height(I))]);
 xlabel('Pitch (nm)','FontSize',10); ylabel('Width (nm)','FontSize',10);
 s1.MarkerEdgeColor = 'k';
s1.MarkerFaceColor = 'r';
set(f1, 'name', 'Per measurement', 'numbertitle', 'off');
set(f1, 'units', 'normalized', 'outerposition', [0.75 0.55 .25 .45]);

f2=figure;
h1=histogram([I{:,1}],100,'FaceColor', 'r');
title(['Morphology per measurement. N=', num2str(height(I))]); 
xlabel('Pitch (nm)','FontSize',10); ylabel('Recurrence','FontSize',10);
set(f2, 'name', 'Per measurement', 'numbertitle', 'off')
set(f2, 'units', 'normalized', 'outerposition', [0.75 0.075 .25 .45]);  
    
% Per fibril
f3=figure;
C = linspace(1,10,length(Fs(:,1)));
s2=scatter(Fs(:,1),Fs(:,3),[],C,'filled');
s2.MarkerEdgeColor = 'k';
title(['Morphology per fibril. N=', num2str(totalFibrils)]);
 xlabel('Pitch (nm)','FontSize',10); ylabel('Width (nm)','FontSize',10);
set(f3, 'name', 'Per fibril', 'numbertitle', 'off')
set(f3, 'units', 'normalized', 'outerposition', [0.5 0.55 .25 .45]);

f4=figure;
h2=histogram(Fs(:,1),50,'FaceColor', 'r');
title(['Morphology per fibril. N=', num2str(totalFibrils)]);
 xlabel('Pitch (nm)','FontSize',10); ylabel('Recurrence (fibril)','FontSize',10);
set(f4, 'name', 'Per fibril', 'numbertitle', 'off')
set(f4, 'units', 'normalized', 'outerposition', [0.5 0.075 .25 .45]);     

% Standard deviation per fibril
f5=figure;
C = linspace(1,10,length(Fs(:,1)));
s3=scatter(Fs(:,1),Fs(:,2),[],C,'filled');
   s3.MarkerEdgeColor = 'k';
   title(['Standard deviation vs. Fibril pitch. N=', num2str(totalFibrils)]);
 xlabel('Pitch (nm)','FontSize',10); ylabel('StdDev (nm)','FontSize',10);
set(f5, 'name', 'Per fibril', 'numbertitle', 'off')
set(f5, 'units', 'normalized', 'outerposition', [0.25 0.55 .25 .45]);  

%% Per fibril and stdev - 3D-scatter plot
f6=figure;
C=cool(height(Fs));
h=scatter3(Fs(:,2),Fs(:,1),Fs(:,3),[],C,'filled','MarkerEdgeColor',[.8 .8 .8], 'SizeData',50 );

title(['Morphology per fibril vs. Std deviation of pitch distance. N=', num2str(totalFibrils)]);
 xlabel('Stdev (nm)','FontSize',10);
 zlabel('Width (nm)', 'FontSize', 10);
 ylabel('Pitch (nm)','FontSize',10);
 set(gca, 'Ydir','reverse', 'Xdir', 'reverse','XminorGrid', 'on', 'YminorGrid', 'on','ZminorGrid', 'on',...
     'MinorGridLineStyle',':','MinorGridColor',[.9 .9 .91], 'MinorGridAlpha',0.25,...
     'GridColor','w','GridAlpha',0.5,'Color', [.15 .15 .15],...
     'LineWidth', 0.5);
 view(-57.3,12);
set(f6, 'name', 'Per fibril', 'numbertitle', 'off')
set(f6, 'units', 'normalized', 'outerposition', [0.25 0.075 .25 .45]);

[x1, y1, z1]=sphere(10);
sizeBall=0.2;
x1=sizeBall*x1; y1=sizeBall*y1; z1=sizeBall*z1;

%% 3D scatter plot with spheres

f7=figure;
set(f7, 'name', 'Per fibril', 'numbertitle', 'off'); set(f7, 'units', 'normalized', 'outerposition', [0 0.55 .35 1]);

subplot(2,1,1);
C=cool(height(Fs));
% h=scatter3(Fs(:,2),Fs(:,1),Fs(:,3),[],C,'filled','MarkerEdgeColor',[.8 .8 .8], 'SizeData',100 );
title(['Morphology per fibril vs. Std deviation of pitch distance. N=', num2str(totalFibrils)]);
 xlabel('Stdev (nm)','FontSize',10); zlabel('Width (nm)', 'FontSize', 10); ylabel('Pitch (nm)','FontSize',10);
 set(gca, 'Ydir','reverse', 'Xdir', 'reverse','XGrid', 'on', 'YGrid', 'on','ZGrid', 'on',...
     'XminorGrid', 'on', 'YminorGrid', 'on','ZminorGrid', 'on',...
     'MinorGridLineStyle',':','MinorGridColor',[.9 .9 .91], 'MinorGridAlpha',0.25,...
     'GridColor',[.8 .8 .8],'GridAlpha',0.5,'Color', [.25 .25 .25],...
     'LineWidth', 0.5);
ax = gca; hold(ax, 'on');

[x1, y1, z1]=sphere(10); sizeBall=0.5; x1=sizeBall*x1; y1=sizeBall*y1; z1=sizeBall*z1;
% Loop through each data element.
for i = 1:height(Fs(:,2))
    if Fs(i,2)==0
        surf(gca, Fs(i,2)+x1, Fs(i,1)+y1, Fs(i,3)+z1, 'FaceColor', 'r', 'EdgeColor', 'none');
    else surf(gca, Fs(i,2)+x1, Fs(i,1)+y1, Fs(i,3)+z1, 'FaceColor', [C(i,:)], 'EdgeColor', 'none'); end
end
h = findall(ax, 'Type', 'surface'); set(h,'FaceLighting','phong', 'AmbientStrength',0.7); light('Position',[1 0 0.5],'Style','infinite');
view(3); set(gca, 'cameraposition',[-55,1990,64],'view',[-16,26],...
    'cameratarget', [10,125,11]);

subplot(2,1,2);
C=cool(height(Fs));
% h=scatter3(Fs(:,2),Fs(:,1),Fs(:,3),[],C,'filled','MarkerEdgeColor',[.8 .8 .8], 'SizeData',100 );
 xlabel('Stdev (nm)','FontSize',10); zlabel('Width (nm)', 'FontSize', 10); ylabel('Pitch (nm)','FontSize',10);
 set(gca, 'Ydir','reverse', 'Xdir', 'reverse','XGrid', 'on', 'YGrid', 'on','ZGrid', 'on',...
     'XminorGrid', 'on', 'YminorGrid', 'on','ZminorGrid', 'on',...
     'MinorGridLineStyle',':','MinorGridColor',[.9 .9 .91], 'MinorGridAlpha',0.25,...
     'GridColor',[.8 .8 .8],'GridAlpha',0.5,'Color', [.6 .6 .6],...
     'LineWidth', 0.5);
ax = gca; hold(ax, 'on');

[x1, y1, z1]=sphere(10); sizeBall=0.5; x1=sizeBall*x1; y1=sizeBall*y1; z1=sizeBall*z1;
% Loop through each data element.
for i = 1:height(Fs(:,2))
    if Fs(i,2)==0
        surf(gca, Fs(i,2)+x1, Fs(i,1)+y1, Fs(i,3)+z1, 'FaceColor', 'r', 'EdgeColor', 'none');
    else surf(gca, Fs(i,2)+x1, Fs(i,1)+y1, Fs(i,3)+z1, 'FaceColor', [C(i,:)], 'EdgeColor', 'none'); end
end
h = findall(ax, 'Type', 'surface'); set(h,'FaceLighting','phong', 'AmbientStrength',0.7); light('Position',[1 0 0.5],'Style','infinite');
view(3);
set(gca, 'cameraposition',[-181,515,48],'view',[-80,15],...
    'cameratarget', [11,121,12]);

%% Final report

fprintf('\n \n >>> A total of %d', totalFibrils); fprintf(' fibrils were measured: \n \n');
fprintf(' Average crossover distance was %f', nanmean(nonzeros(Fs(:,1))))
fprintf(' +/- %f', nanmean(nonzeros(Fs(:,2)))); fprintf(' nanometers. \n');
fprintf(' Average width was %f', nanmean(nonzeros(Fs(:,3))))
fprintf(' +/- %f', nanmean(nonzeros(Fs(:,4)))); fprintf(' nanometers. \n');
fprintf(' Assuming C2 symmetry with a rise of 4.75 A, average twist is +/- %f', (0.475*180)/nanmean(nonzeros(Fs(:,1)))); 
fprintf(' degrees \n');
fprintf(' Assuming 2n-screw symmetry with a rise of 2.375 A, average twist is +/- %f', (360-(0.475*180)/nanmean(nonzeros(Fs(:,1))))/2); 
fprintf(' degrees \n \n');

% Looking only at fibrils below 100 nm pitch
counter=1; tot=0;
for a=1:height(Fs)
    if Fs(a,1)<=100 
        Fshort(counter,:)=Fs(a,:);
        counter=counter+1;end
end

fprintf('\n >>> Of all fibrils, %d had a cross=over distance below 100 nm \n', height(Fshort));
fprintf('\n Average crossover distance for those was %f', nanmean(nonzeros(Fshort(:,1))));
fprintf('  +/-  %d nanometers \n', nanmean(nonzeros(Fshort(:,2))));
fprintf(' Average width was %d   +/-   %d nanometers. \n', nanmean(nonzeros(Fshort(:,3))), nanmean(nonzeros(Fshort(:,4))));
fprintf(' Assuming C2 symmetry with 4.75 A rise, average twist is +/- %d degrees \n', (0.475*180)/nanmean(nonzeros(Fshort(:,1))));
fprintf(' Assuming 2n-screw symmetry with 2.375 A rise, average twist is +/- %d degrees \n',(360-(0.475*180)/nanmean(nonzeros(Fshort(:,1))))/2);






 