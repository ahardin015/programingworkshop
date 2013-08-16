clear all %clean out the Workspace

%load in necessary matrix data for plotting  
load('lat.mat');
load('lon.mat');
load('lat_ref.mat');
load('lon_ref.mat');
load('nwsmap.mat'); %for plotting
%Set interpolated grid lat and lon
newlat=lat_ref;
newlon=lon_ref;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defined contstant Values=================================================
original_grid=3/7; %originally 3/7 for model euclidean, but set to 0 for lat/lon
interpolated_grid=9/7;
sigma=1.25; %smoothing factor into gaussian smoothing (higher smoothes more)
masking_threshold=25; %Note, this should be lower than ci_threshold
ci_threshold=35;
object_threshold=2; %Number of grid points to search(of interpolated grid spacing)
match_spatial_thresh=10; %Number of grid space distance units (depedent on interp grid)
match_temporal_thresh=2; %Number of time step units (5 min steps)
num_match_thresh=5; %Number of times object needs to be matched in time to be flagged true CI object
sim_time=1; %12 frames/hr 
time='2200'; %Starting Z time string (must be 4 digits)
date='0621'; %date string for saving vars at end
%==========================================================================

%open simulated dBz file===================================================
dbzname=['dbz' date time '_mod'];
load(dbzname)
%==========================================================================

%Allocation================================================================
master_obslist=NaN(1,19);
object_grid=NaN(size(newlat,1),size(newlat,2),sim_time);
%==========================================================================
% yprevstart=1; %these variables will keep track of the index of object list from previous time step
yprevend=1;

dbzall=eval(dbzname);

%Time Dimension
for t=1 %modified to read in t=37 only (22z)
  
%Load in model dBz field at time t
dbz=dbzall(:,:,t);
% dbz=dbz05292000_mod(:,:,t);
    
%Take Out negative dBz values=============================================
dbz_new=zeros(size(dbz,1),size(dbz,2));
for i=1:size(dbz,1)
    for j=1:size(dbz,2)
        if isfinite(dbz(i,j))==1 && dbz(i,j) >= 0
            dbz_new(i,j)=dbz(i,j);
        end
    end
end

%Figure 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
contourf(dbz_new)
shading flat
colorbar
title(['Original simulated -10 C dBz field at ' time 'Z'])
hold on
colormap(nwsmap)
caxis([0 75])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Up-interpolate to a coarser grid ...')
%Call fit_to_gridb_test====================================================
%Put 3/7km eucl. grid to coarser eucl. grid for computational efficiency
[newdat] = fit_to_gridb_test(double(dbz_new),double(lat),double(lon),lat_ref,lon_ref,original_grid,interpolated_grid,true);
%==========================================================================

%Figure 2%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
contourf(newlon,newlat,newdat)
shading flat
colorbar
title(['Up-interpolated simulated -10 C dBz field at ' time 'Z'])
hold on
colormap(nwsmap)
caxis([0 75])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Smoothing the data ...')
%Call gaussian smooth (Convolution)=======================================
% [pdf] = gaussian_smooth(mat,fact)
[pdf] = gaussian_smooth(newdat,sigma);

%Figure 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
contourf(newlon,newlat,pdf)
shading flat
hold on
colorbar
title(['Gaussian smoothed simulated -10 C dBz field at ' time 'Z'])
hold on
colormap(nwsmap)
caxis([0 75])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
disp('Masking out values that do not meet CI criteria ...')
%Mask out values =========================================================
dbz_smooth=pdf;
dbz_masked=NaN(size(dbz_smooth,1),size(dbz_smooth,2));
for i=1:size(dbz_smooth,2)
    for j=1:size(dbz_smooth,1)
        if dbz_smooth(j,i)>=masking_threshold && newdat(j,i)>=ci_threshold
%Insert original dBz values back in=======================================            
            dbz_masked(j,i)=newdat(j,i);
        end
    end
end

%Figure 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
contourf(newlon,newlat,dbz_masked)
shading flat
hold on
colorbar
title(['Masked Gaussian smoothed simulated -10 C dBz field at ' time 'Z'])
hold on
colormap(nwsmap)
caxis([0 75])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
disp('group and define CI objects ...')
%Group into Objects========================================================
dbz_obj=dbz_masked;

x=1:1:size(dbz_obj,1);
y=1:1:size(dbz_obj,2);

[ogridy,ogridx]=meshgrid(y,x);
dbz_res=reshape(dbz_obj,size(dbz_obj,1)*size(dbz_obj,2),1);
% dbz_inds=find(~isnan(dbz_obj));
dbz_inds=find(dbz_obj>0);
dbz_list=dbz_res(dbz_inds);
%reshape x and y to find x,y indices
x_res=reshape(ogridx,size(ogridx,1)*size(ogridx,2),1);
y_res=reshape(ogridy,size(ogridy,1)*size(ogridy,2),1);

x_point=x_res(dbz_inds);
x_list=x_point(dbz_list>0); %Finds non-zero entries of dbz_list

y_point=y_res(dbz_inds);
y_list=y_point(dbz_list>0);

%Define Objects============================================================
%==========================================================================

%beginning index in master_obslist of last time step was clsuter_ctin
%before it is reset below


%Call subroutine objects
if ~isempty(x_list)%if not empty
    
    if exist('cluster_ctin','var')==0 %if does not exist
        cluster_ctin=0; %declare it
    else
        cluster_ctin=cluster_ctout; %else set it to last count
    end
    
cluster_ctin

[observed_object_list, cluster_ctout] = objects(x_list,y_list,cluster_ctin,object_threshold,interpolated_grid);

cluster_ctout

%returns a list of corresponding object numbers to x_list and y_list
%convective reflectivity points

%send dbz,x, and y _lists to shape orientation
% [slope, xlength, ylength] = shape_orientation(x_list,y_list,dbz_list);
    %slope: angle of semi-major axis relative to x-axis
    %xlength: length of ellipse semi-major axis
    %ylength: length of ellipse semi-minor axis

% object_grid(:,:,t)=NaN(size(newlon,1),size(newlat,2));
%create a 2D matrix of object number for convectively active grid points
for r=1:length(observed_object_list)
    if isnan(observed_object_list(r))==0
    object_grid(x_list(r),y_list(r),t)=observed_object_list(r);
    %make sure zeros from objects that didn't reach size criteria are set
    %back to NaN
    end
end

% object_grid_test=zeros(size(object_grid,1),size(object_grid,2),sim_time);
% object_grid_test(:,:,t)=object_grid(:,:,t);
% bad_inds=find(isnan(object_grid_test));
% object_grid_test(bad_inds)=0;
% clear bad_inds

%==========================================================================

%Matching Objects==========================================================

%First set of objects====================================================
%Add first group of objects into master list with attributes calculted,
%then later objects will either be added in with attributes or if matched
%with previous object, tally into temporal CI threshold
 if cluster_ctin==0 && cluster_ctout>cluster_ctin 
    
    for ob=1:cluster_ctout
    
    [rows, cols, refs]=find(object_grid(:,:,t)==ob);
    master_obslist(ob,1)=ob; %Object number
    master_obslist(ob,2)=0; %CI flag (becomes 1 when matched 3 times)
    master_obslist(ob,3)=0; %Match count in time (tally)
    master_obslist(ob,4)=str2double(time); %First time (Z time)
    master_obslist(ob,5)=str2double(time); %Latest time (Z time)
    %Area average object center
    master_obslist(ob,6)=sum(cols)/length(cols); %cols_center(ob)
    master_obslist(ob,7)=sum(rows)/length(rows); %rows_center(ob)
    master_obslist(ob,8)=(length(rows))*(interpolated_grid^2); %object_area(ob)
    %add in time t attribute (9/21/12)
    master_obslist(ob,9)=t; %latest time (integer count)
    %calc lat/lon center from new array from new_object
    new=master_obslist; %trick for not replacing all "new" arrays with "master_oblist" arrays
    master_obslist(ob,10)=lon_ref(floor(new(ob,7)),floor(new(ob,6))) + ((new(ob,6)-floor(new(ob,6)))*(lon_ref(floor(new(ob,7)),floor(new(ob,6))+1)-lon_ref(floor(new(ob,7)),floor(new(ob,6)))));%lon_center (from cols_center)
    master_obslist(ob,11)=lat_ref(floor(new(ob,7)),floor(new(ob,6))) + ((new(ob,7)-floor(new(ob,7)))*(lat_ref(floor(new(ob,7))+1,floor(new(ob,6)))-lat_ref(floor(new(ob,7)),floor(new(ob,6)))));%lat_center (from rows_center)
    master_obslist(ob,12)=0; %time of CI declaration, not yet occurred here
    master_obslist(ob,13)=sum(cols)/length(cols); %sum of cols_center through CI declaration 
    master_obslist(ob,14)=sum(rows)/length(rows); %sum of rows_center through CI declaration
    master_obslist(ob,15)=t; %first time object identified
    master_obslist(ob,16)=sum(cols)/length(cols); %first CA time cols loc
    master_obslist(ob,17)=sum(rows)/length(rows); %first CA time rows loc
%Weighted Object Center
%     attributes(ob,4,t)=sum(rows.*refs)/sum(refs); %rows_center_weighted(ob)
%     attributes(ob,5,t)=sum(cols.*refs)/sum(refs);
%     %cols_center_weighted(ob)    
clear rows cols refs new
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear latlist lonlist numlist
%======================================================================== 
%Handle New Object======================================================
 elseif cluster_ctin>0 && cluster_ctout>cluster_ctin
        %ending index of master_obslist from last time step is length of
        %array before adding new objects at current time step
   closest_match=zeros(cluster_ctout-cluster_ctin,2);
   
    for ob=cluster_ctin+1:cluster_ctout %go through objects of time t
    [rows, cols, refs]=find(object_grid(:,:,t)==ob);
    %Temporary attributes
    cols_center=sum(cols)/length(cols);
    rows_center=sum(rows)/length(rows);
        
    %Spatial distances that are within temporal threshold boundaries are
    %calculated
    time_inds=find(( (abs(master_obslist(:,9)-t)<=match_temporal_thresh) & (master_obslist(:,9)-t)~=0));
        if size(time_inds,1)>0    
    
    colsdist=master_obslist(time_inds,6) - cols_center;
    rowsdist=master_obslist(time_inds,7) - rows_center;
    dist(:,1)=time_inds;
    dist(:,2)=sqrt(colsdist(:).^2 + rowsdist(:).^2);
    
    %Temporal distance (absolute diff since using a less than statement)
%    timediff=abs(master_obslist(yprevstart:yprevend,9)-t);

    row_inds=find(((dist(:,2)<=match_spatial_thresh)));%&(timediff<=match_temporal_thresh));
    
    %handle if 1 or more close matches
    if length(row_inds) > 1
        small=min(dist(row_inds));
        clear row_inds
        row_inds=find(dist==small);
        closest_match(ob-cluster_ctin,1)=dist(row_inds,1); %linear indice of closet match of past object
        closest_match(ob-cluster_ctin,2)=dist(row_inds,2); %dist from closet match of past object  
    elseif length(row_inds)==1 
     %write in closest matched indices and distance
    closest_match(ob-cluster_ctin,1)=dist(row_inds,1); %linear indice of closet match of past object
    closest_match(ob-cluster_ctin,2)=dist(row_inds,2); %dist from closet match of past object  
    end
        end
    clear row_inds rows cols refs dist  
    end
%========================================================================

    %Look for multiple matches to a single past object
    %(i.e at t=2 3 objects closely match an object from time t=1)
    for kk=1:size(closest_match,1)
       if closest_match(kk,1)~=0
        diff=closest_match(:,1)-closest_match(kk,1); %difference of kk object inds match to all others
        inds=find(diff==0); %indices of other objects that have same match in past objects but not already set to zero
        
        if length(inds) > 1 %a match besides itself to a past object
        competitors(1,1)=closest_match(kk,2); %current ob num kk first in list
        competitors(2:length(inds),1)=closest_match(inds(2:length(inds)),2); %others obs added
        closest=find(closest_match(:,2)==min(competitors)); %find who is closest of list 
        for tt=1:length(inds) %go through length of inds
            if inds(tt) ~= closest %if not closet, zero out its match from past object
                closest_match(inds(tt),1:2)=0; %check this and now add in another call to new_objects, or just move it down
            end
        end  
    clear diff ind competitors closest
        end
       end
    end
%match closest objects to past
        for ob=cluster_ctin+1:cluster_ctout
            
            [rows, cols, refs]=find(object_grid(:,:,t)==ob);
             
            if closest_match(ob-cluster_ctin,1)==0 %new object

            clear new
            [new] = new_object(rows,cols,refs,interpolated_grid,time,ob,t);    
            newloc=size(master_obslist,1)+1;
            master_obslist(newloc,1:9)=new;
            master_obslist(newloc,1)=newloc;
            %calc lat/lon center from new array from new_object
            master_obslist(newloc,10)=lon_ref(floor(new(1,7)),floor(new(1,6))) + ((new(1,6)-floor(new(1,6)))*(lon_ref(floor(new(1,7)),floor(new(1,6))+1)-lon_ref(floor(new(1,7)),floor(new(1,6)))));%lon_center (from cols_center)
            master_obslist(newloc,11)=lat_ref(floor(new(1,7)),floor(new(1,6))) + ((new(1,7)-floor(new(1,7)))*(lat_ref(floor(new(1,7))+1,floor(new(1,6)))-lat_ref(floor(new(1,7)),floor(new(1,6)))));%lat_center (from rows_center)
            master_obslist(newloc,12)=0; %not a CI object yet
            master_obslist(newloc,13)=new(1,6); %sum of cols_center through CI declaration 
            master_obslist(newloc,14)=new(1,7); %sum of rows_center through CI declaration
            master_obslist(newloc,15)=t; %first time object identified
            master_obslist(newloc,16)=new(1,6); %first CA time cols loc
            master_obslist(newloc,17)=new(1,7); %first CA time rows loc
            
            else %match to closest object
                
            row_inds=closest_match(ob-cluster_ctin,1);
            
            %CASE1: whether or not flagged CI, update latest attributes    
            master_obslist(row_inds,3)=master_obslist(row_inds,3)+1; %update count tally
            master_obslist(row_inds,5)=str2double(time); %z time string
            master_obslist(row_inds,6)=sum(cols)/length(cols); %latest object centroid col
            master_obslist(row_inds,7)=sum(rows)/length(rows); %latest object centroid row
            master_obslist(row_inds,9)=t; %latest time ,t, matched for finding objects within temporal range
           
            %CASE2: object not yet flagged CI, accumulate col,row center and area prior to CI 
            if master_obslist(row_inds,3)<=num_match_thresh     
                master_obslist(row_inds,8)=master_obslist(row_inds,8)+(length(rows))*(interpolated_grid^2); %sum area
                master_obslist(row_inds,13)=master_obslist(row_inds,13)+(sum(cols)/length(cols)); %sum centroid col
                master_obslist(row_inds,14)=master_obslist(row_inds,14)+(sum(rows)/length(rows)); %sum centroid row  
            end
            
            %CASE3: object now a CI object, calculate average CI attributes
            if master_obslist(row_inds,3)==num_match_thresh
                master_obslist(row_inds,2)=1; %Flagged true CI object
                master_obslist(row_inds,8)=master_obslist(row_inds,8)/(num_match_thresh+1); %avg accum CI object area
                %add in lat/lon centers for confirmed CI objects (9/27/12)
                master_obslist(row_inds,10)=lon_ref(floor(master_obslist(row_inds,7)),floor(master_obslist(row_inds,6))) + ((master_obslist(row_inds,6)-floor(master_obslist(row_inds,6)))*(lon_ref(floor(master_obslist(row_inds,7)),floor(master_obslist(row_inds,6))+1)-lon_ref(floor(master_obslist(row_inds,7)),floor(master_obslist(row_inds,6)))));%lon_center (from cols_center)
                master_obslist(row_inds,11)=lat_ref(floor(master_obslist(row_inds,7)),floor(master_obslist(row_inds,6))) + ((master_obslist(row_inds,7)-floor(master_obslist(row_inds,7)))*(lat_ref(floor(master_obslist(row_inds,7))+1,floor(master_obslist(row_inds,6)))-lat_ref(floor(master_obslist(row_inds,7)),floor(master_obslist(row_inds,6)))));%lat_center (from rows_center)
                master_obslist(row_inds,12)=t; %Time of CI declaration
                master_obslist(row_inds,13)=master_obslist(row_inds,13)/(num_match_thresh+1); %avg accum centroid col
                master_obslist(row_inds,14)=master_obslist(row_inds,14)/(num_match_thresh+1); %avg accum centroid row
                master_obslist(row_inds,18)=master_obslist(row_inds,6); %cols loc at CI declaration
                master_obslist(row_inds,19)=master_obslist(row_inds,7); %rows loc at CI declaration
            end
            end
            clear rows cols refs
        end
        clear closest_match

% if strcmp(time(1,3:4),'30')==1 || strcmp(time(1,3:4),'00')
%only print objects if >1 defined object for contourfm
%     if (last-yprevend) > 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% end  
clear latlist lonlist numlist    
 end
%==========================================================================    
end %if xlist not empty (if no xlist, no objects, skip to next time)
last=size(master_obslist,1);

%Figure 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
latlim=[39.101150512695310 44.463722229003906];
lonlim=[-105.7975921630859 -99.627136230468750];
ax = usamap(latlim, lonlim);
latlim = getm(ax, 'MapLatLimit');
lonlim = getm(ax, 'MapLonLimit');
states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
contourfm(newlat,newlon,object_grid(:,:,t))
geoshow(ax, states, 'FaceColor', 'none')
hold on
latlist=master_obslist(yprevend+1:last,11);
lonlist=master_obslist(yprevend+1:last,10);
numlist=num2cell(master_obslist(yprevend+1:last,1));
plotm(latlist(:),lonlist(:),'k+')
hold on
textm(latlist(:),lonlist(:),numlist(:),'VerticalAlignment','bottom','HorizontalAlignment','right')
hold on
title(['Simulated objects at ' time 'Z (+ = new object; else matched)'])
colorbar
% print(gcf, '-dpdf', [date time 'modobj.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     end
%Figure 6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
latlim=[39.291255950927734 44.463722229003906];
lonlim=[-105.7975921630859 -99.627136230468750];
ax = usamap(latlim, lonlim);
latlim = getm(ax, 'MapLatLimit');
lonlim = getm(ax, 'MapLonLimit');
states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
geoshow(ax, states, 'FaceColor', 'none')
contourfm(newlat,newlon,newdat)
hold on
shading flat
title(['Up-interpolated simulated -10 C dBz field at ' time 'Z'])
hold on
colormap(nwsmap)
colorbar
caxis([0 75])
geoshow(ax, states, 'FaceColor', 'none') %run geoshow again
% print(gcf, '-dpdf', [date time 'moddbz.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dbz
clear dbz_list
clear x_list
clear y_list
clear observed_object_list

    %advance time counter 5 minutes accounting for hour changes
%     mm=time(3:4);
%     hh=time(1:2);
%     
%     mm=str2double(mm);
%     mm=mm+5;
%     if mm == 60
%        mm=0;
%        hh=str2double(hh);
%        hh=hh+1;
%        hh=num2str(hh);
%        if length(hh) < 2
%            hh = ['0' hh];
%        end
%     end
%     mm=num2str(mm);
%     if length(mm) < 2
%      mm = ['0' mm];
%     end
%     %new time string
%     time = [hh mm];
%     
%     yprevend=size(master_obslist,1);
end
