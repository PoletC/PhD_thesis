%% Cod for Polet %%
clear all
close all
disp(['Now analysing:',pwd])

tb_analysed = dir('*.tif');
minimum_cell_pix_size = 1000; %usually size of cells is above 1K pixels in images;
Sensitivity= 0.3; %Vary this parameter if you think the fit is poor; usually its not larger than 0.5!
mkdir(pwd,'AnalysedData'); %create a folder to save all the output in
Disk1Size = 5;%change this to make the sectioningmore/less. Very senstive!!!
Disk2Size = 15;%change this to  fill the outlines more/less aggressive
% Disk3Size = 10;%change this to fill the outlines more/less (not used
% anymore!)
SaveImages = true; %set this to true (SaveImages = true;) to save shown figures

set(0,'DefaultFigureVisible','off'); %comment this to see the figures!

for i=1:length(tb_analysed)
    
    file = tb_analysed(i).name;
    path = tb_analysed(i).folder;
    %we first locate the multipage tiff we want to analyse and load it into the
    %matlab memory
    img_info = imfinfo(file);
    tiff_sizex = img_info.Width;
    param.width=tiff_sizex-1;%save this for future reference, depends on binning!
    tiff_sizey = img_info.Height;
    param.height=tiff_sizey-1;%save this for future reference, depends on binning!
    tiff_sizez = size(imfinfo(file),1);

     
    %savename creation
    save_name_string = strsplit(file,'.');
    save_name_string = fullfile('AnalysedData',save_name_string{1}); %the bit before the point at the end
    
    %load in the image
    I = imread(fullfile(path, file));

    bw = imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',Sensitivity);
%     bw = imbinarize(I,'global'); binarize image 

    imshow(bw)
    pause(0.5)
    
    bw2 = bw; %invert if ForegroundPolarity is 'bright'
    
    %first we fill the voids
    se = strel('disk',Disk1Size);
    bw_trial1 = imclose(bw2,se);
    imshow(bw_trial1);
    pause(0.5)
    
    %then we remove objects smaller than 10 pixels in radius
    se2 = strel('disk',Disk2Size);
    bw_trial3 = imopen(bw_trial1,se2);
    imshow(bw_trial3)
    pause(0.5)
    
    %then we fill the remaining voids
%     se = strel('disk',Disk3Size);
%     bw_trial4 = imclose(bw_trial3,se);
%     imshow(bw_trial4);
%     pause(0.5)    

%     final_analysis = ~bw_trial4;
    final_analysis = ~bw_trial3;

    %show the final masked image and overlay the found areas using
    %regionprops
    imshow(final_analysis)

    stats = regionprops('table',final_analysis,'Centroid',...
        'MajorAxisLength','MinorAxisLength','Area','Orientation','Perimeter');
    centers = stats.Centroid;
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = diameters/2;
    hold on
    viscircles(centers,radii);
    hold off

    pause(2)
    close all

    all_elements = round(table2array(stats));
    %remove smallest, area larger than 1000 pixels at least
    final_elements = all_elements(all_elements(:,1)>minimum_cell_pix_size,:);
    
    imshow(final_analysis)
    centers2 = final_elements(:,2:3);
    diameters2 = mean(final_elements(:,4:5),2);
    radii2 = diameters2/2;
    hold on
    if length(centers2(:,1))>0
        for j=1:length(centers2(:,1))
            h = drawellipse('Center',centers2(j,:),'SemiAxes',final_elements(j,4:5)./2,'RotationAngle',final_elements(j,6));
        end
    end
    %viscircles(centers2,radii2);
    hold off
    
    if SaveImages==true
        print([save_name_string,'Selected_Cells_For_Analysis.tif'],'-dtiff')
    end
    
    pause (1)
    close all

 %   imshowpair(I,final_analysis,'montage')
    imshowpair(I,final_analysis)
    title('Original image + analysed areas of cells')

    if SaveImages==true
        print([save_name_string,'after_analysis.tif'],'-dtiff')
        axis on
        savefig([save_name_string,'after_analysis.fig'])
    end
    close all

    imshow(I,[min(double(I(:))) max(double(I(:)))])
    if length(centers2(:,1))>0
        for j=1:length(centers2(:,1))
            h = drawellipse('Center',centers2(j,:),'SemiAxes',final_elements(j,4:5)./2,'RotationAngle',final_elements(j,6));
        end
    end
    if SaveImages==true
        print([save_name_string,'after_analysis_on_original.tif'],'-dtiff')
        savefig([save_name_string,'after_analysis_on_original.fig'])    
    end
    
    final_elements_forsaving = num2cell(final_elements);
    final_elements_forsaving2 = vertcat({'All units are pixels of the numbers below!','','','','','',''},{'Area','Centroidx','Centroidy','MajorAxisLength','MinorAxisLength','Orientation','Perimeter'},final_elements_forsaving);
    writecell(final_elements_forsaving2,[save_name_string,'_dataCells.xlsx'])
    
    %save workspace
    save([save_name_string,'_workspace.mat'])

    close all
    disp(num2str(i))
    clearvars -except tb_analysed minimum_cell_pix_size Sensitivity Disk1Size Disk2Size Disk3Size SaveImages

end

set(0,'DefaultFigureVisible','on');

