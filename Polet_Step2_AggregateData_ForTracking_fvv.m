%% Code for Polet %%

% First, navigate the Working Directory (pwd) to the AnalysedData folder
% made in step 1!!!
clear all
close all
disp(['Now analysing the .mat files in :',pwd])

% set this to false to create the large image files!
%Memsmall = false; % in case you have more than 20 GB RAM...
Memsmall = true;

% find all the workspaces
tb_analysed_mats = dir('*workspace.mat');

%ensure that plots are plotted
set(0,'DefaultFigureVisible','on'); %comment this to see the figures!

%Load in the first step: then use this info to praollocate space in memory
file = tb_analysed_mats.name;
TBloadedObject = matfile(file);
varlist = who(TBloadedObject);
% load in I,final_analysis, final_elements
load(file,'I','final_analysis','final_elements','minimum_cell_pix_size');
disp(['The size cutoff used in this dataset is: ',num2str(minimum_cell_pix_size),' pixels.'])

% determine how many images (timesteps) are in this folder
numsteps = length(tb_analysed_mats);

%preallocate space for the xyzt list (for tracking) and the images stack
XYT_list_for_tracking = [];

if Memsmall == false %
    I_all_timesteps = nan([size(I) numsteps]);
    FinalMask_all_timesteps = nan([size(I) numsteps]);
end


Final_Objects_all_timesteps = cell(numsteps,1);
Stats_AllObjects_all_timesteps = cell(numsteps,1);

for k=1:numsteps
    
    file = tb_analysed_mats(k).name;
    path = tb_analysed_mats(k).folder;
    
    TBloadedObject = matfile(file);
    varlist = who(TBloadedObject);
    %load in the analysed data (only the positions and image at this
    %timestep
    load(file,'I','final_analysis','final_elements','stats');
    
    %save the stats in a struct for each timestep
    Final_Objects_all_timesteps(k,1) = {final_elements};
    Stats_AllObjects_all_timesteps(k,1) = {stats};
    
    %assemble an image stack of all the timesteps (this is the slow part,
    %feel free to cancel if you want. Very memory intensive!!!
    if Memsmall == false %
        I_all_timesteps(:,:,k) = I;
        FinalMask_all_timesteps(:,:,k) = final_analysis;
    end
    
    %assemble a list of the position in the format xyt for the tracking
    %code
    
%   final_elements numer order: 'Area','Centroidx','Centroidy',...
%   'MajorAxisLength','MinorAxisLength','Orientation','Perimeter'
    positions_thisstep = final_elements(:,2:3); %x, y position
    positions_thisstep_plustime = [positions_thisstep zeros(length(positions_thisstep),1)+k];
    %add the timepoint onto the array
    
    XYT_list_for_tracking = [XYT_list_for_tracking;positions_thisstep_plustime];
    
    disp(['Now loading timestep: ',num2str(k)])
end

%% run the tracking code %%

maxDisplacement = 300;%max displacement, should be smaller then average interparticle distance
% I assumed a 300 pixels (which is 2x the cell size approximately) would be a
% good estimate
param.mem = 20; %steps in a row the particle can be missed but still be part of a track :)
param.dim = 2; %for the two dimensional data(X,Y)
param.good = ceil(numsteps/10); %the minimum track length, 10% of total track
param.quiet = 0; %set to 1 to suppress text

tracks = track(XYT_list_for_tracking,maxDisplacement,param);
%tracks is x,y,timestep,CellID
tracks_SortTime = sortrows(tracks,3);
% now we want to add the other regionprops data to the found tracks

% make the thing we save all the data in:
tracks_complete = tracks;

Colorscale_steps = max(tracks(:,4));%the maximum beadID found
colormap_forAggr = cbrewer2('div','Spectral',Colorscale_steps,'spline');
colormap_forAggr = colormap_forAggr./max(colormap_forAggr);

mkdir('Movietracks')

if Memsmall == false %
    MeanMax = mean(max(max(I_all_timesteps(:,:,:))));
    %compute a max value for intensity to plot later. In case of small
    %memory, the reference image is used so this is not needed.
end

for l=1:numsteps
    
    RegPropsData = Final_Objects_all_timesteps{l,1};%load back in the data
    TracksData = tracks_SortTime(tracks_SortTime(:,3)==l,:); % the timestep data from the track
    
    for p=1:max(size(TracksData(:,1)))
        % find the other info (Cell size, etc in the RegPropsData)
        % we use the x coordinate
        rowRegProp = find(TracksData(p,1)==RegPropsData(:,2) & TracksData(p,2)==RegPropsData(:,3));
        % use the found index to get the cellsize etc data back)
        cellinfo_this_timestep = RegPropsData(rowRegProp,:);
        
        %find the original spot back in the tracks
        rowTracksOriginal = find(TracksData(p,1)==tracks(:,1) &TracksData(p,2)==tracks(:,2) &TracksData(p,3)==tracks(:,3));
        
        tracks_complete(rowTracksOriginal,5) = cellinfo_this_timestep(1);% the size of the cell
        tracks_complete(rowTracksOriginal,6) = cellinfo_this_timestep(4);% the ellipse parameter 1
        tracks_complete(rowTracksOriginal,7) = cellinfo_this_timestep(5);% the ellipse parameter 2
        tracks_complete(rowTracksOriginal,8) = cellinfo_this_timestep(6);% the angle of the fit ellipse
        tracks_complete(rowTracksOriginal,9) = cellinfo_this_timestep(7);% the amount of pixels on the border
        
    end
    
     TracksData_complete = tracks_complete(tracks_complete(:,3)==l,:); % the timestep data from the track
   
    
    %savename creation
    save_name_movie_string = fullfile(path,'\Movietracks'); %the bit before the point at the end
    
    % make a plot to show the time evolution!
    figure('units','normalized','outerposition',[0 0 1 1]) % this is to
%     make the image larger on screen
    if Memsmall == false %
        imagesc(I_all_timesteps(:,:,l))
        axis equal
        axis on
        colormap(gray)
        caxis([0 MeanMax/1.2])
    else
        imagesc(I)
        axis equal
        axis on
        colormap(gray)
        clim([0 max(max(I))/1.2])
    end
    
    centers2 = TracksData_complete(:,1:2);
    diameters2 = mean(TracksData_complete(:,6:7),2);
    radii2 = diameters2/2;
    hold on
    if length(centers2(:,1))>0
        for j=1:length(centers2(:,1))
            h = drawcircle('Center',centers2(j,:),'Radius',diameters2(j,1)/2,'Color',colormap_forAggr(TracksData_complete(j,4),:));
%             h = drawellipse('Center',centers2(j,:),'SemiAxes',TracksData_complete(j,6:7)./2,'RotationAngle',TracksData_complete(j,8),'Color',colormap_forAggr(TracksData_complete(j,4),:));
        end
    end
    exportgraphics(gcf,[save_name_movie_string,'\Timestep',sprintf('%04d',l),'.tif'],'Resolution',150)
%     print([save_name_movie_string,'\Timestep',num2str(l),'.tif'],'-dtiff','-r300')
%     savefig([save_name_movie_string,'Timestep',num2str(l),'.fig'])    
    
    close all

    disp(['Completed timestep: ',num2str(l)])
end


%% make plots of the size for each CellID: size and relative size %%

% Colorscale_steps contains the max number of beadIDs (used before)

mkdir('CellsPlotsSizes')
%savename creation
save_name_plots_string = fullfile(path,'\CellsPlotsSizes'); %the bit before the point at the end

for m=1:Colorscale_steps
    
    rows_thisCell = find(tracks_complete(:,4)==m);
    
    timesteps_thisCell = tracks_complete(rows_thisCell,3);
    sizes_thisCell = tracks_complete(rows_thisCell,5);
    %take first 5 values as the mean to make the size increase relative:
    relative_sizes_thisCell = sizes_thisCell./mean(sizes_thisCell(1:5));
    
    Cell_timesteps{m} = {timesteps_thisCell};
    Cell_sizes{m} = {sizes_thisCell};
    Cell_sizes_relative{m} = {relative_sizes_thisCell};
    
    %
    scatter(timesteps_thisCell,sizes_thisCell,'MarkerEdgeColor','k','MarkerFaceColor',colormap_forAggr(m,:))
    box on
    ylim([0 max(tracks_complete(:,5))]);
    xlim([0 max(tracks_complete(:,3))]);
    title(['CellID = ',num2str(m)])
    ylabel(['Size (pixels)'])
    xlabel(['Timestep'])
    exportgraphics(gcf,[save_name_plots_string,'\Cell',num2str(m),'SizeVsTime.tif'],'Resolution',150)
    close all
    
    scatter(timesteps_thisCell,relative_sizes_thisCell,'MarkerEdgeColor','k','MarkerFaceColor',colormap_forAggr(m,:))
    box on
    xlim([0 max(tracks_complete(:,3))]);
    title(['CellID = ',num2str(m),' Relative size'])
    ylabel(['Size (a.u., relative to first 5 obs)'])
    xlabel(['Timestep'])
    exportgraphics(gcf,[save_name_plots_string,'\Cell',num2str(m),'RelSizeVsTime.tif'],'Resolution',150)
    legend
    close all
end

figure;
for m=1:Colorscale_steps
    scatter(cell2mat(Cell_timesteps{m}),cell2mat(Cell_sizes_relative{m}),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceColor',colormap_forAggr(m,:))
    box on
    xlim([0 max(tracks_complete(:,3))]);
    ylabel(['Size (a.u., relative to first 5 obs)'])
    xlabel(['Timestep'])
    legend
    hold off
end
exportgraphics(gcf,[save_name_plots_string,'\AllPlotsRelsizeVsTime.tif'],'Resolution',150)
close all

figure;
for m=1:Colorscale_steps
    scatter(cell2mat(Cell_timesteps{m}),cell2mat(Cell_sizes{m}),'MarkerEdgeColor','k','MarkerEdgeAlpha',0.5,'MarkerFaceColor',colormap_forAggr(m,:))
    box on
    xlim([0 max(tracks_complete(:,3))]);
    ylabel(['Size (a.u., relative to first 5 obs)'])
    xlabel(['Timestep'])
    legend

    hold on
end
exportgraphics(gcf,[save_name_plots_string,'\AllPlotsSizeVsTime.tif'],'Resolution',150)
close all

%% save the workspace %%

if Memsmall == false %
    save_name_alldata_string = fullfile(path,'\AllDataAggrLarge.mat'); %the bit before the point at the end
    save(save_name_alldata_string,'-v7.3')
end
save_name_alldata_string = fullfile(path,'\AllDataAggrSmall.mat'); %the bit before the point at the end
save(save_name_alldata_string,'-v7.3')
%writecell(save_name_alldata_string, fullfile(path,'\AllDataAggrSmall.xlsx'))

%% print matrix
col_header={'Xpos','Ypos','Timestep','CellID','Area','MajorAxisLength','MinorAxisLength','Orientation','Perimeter'};     %Row cell array (for column labels)
output_matrix=[col_header;  num2cell(tracks_complete)];
writecell(output_matrix,[path,'\tracks_complete.xlsx'])

% writematrix("tmp.xlsx", tmp.tracks_complete)

