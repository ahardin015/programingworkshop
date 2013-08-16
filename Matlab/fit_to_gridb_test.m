%This properly trims to original domain box rather than first fit_to_grib
function [newdat] = fit_to_gridb_test(data,lat,lon,lat_ref,lon_ref,dx_old,dx,isuclid)
%input arguments:
% data=double(dbz);
% lat=double(lat);
% lon=double(lon);
% lat_ref=39.2913;
% lon_ref=-105.7976;
% dx_old=3/7;
% dx=5;
% isuclid=true;
%
%Data: original data
%
%lat/lon: original meshgridded (for uclidian grids) or non-meshgridded 
%(for latitude/longitude grids) latitudes and longitudes of data
%
%lat_ref/lon_ref: latitude/longitude of the lower left corner of the grid
%to which the data will be interpolated (for latitude/longitude grids).
%These inputs are only valid when isuclid is set to "false", indicating
%that the grid is a latitude/longitude grid.
%
%dx_old: for uclidian grids, the original grid spacing of the data.  For
%non-uclidain grids, this can be set to zero
%
%dx: grid spacing for grid to which the data is being interpolated
%
%isuclid: is the grid uclidian (true) or latitude/longitude(false)

%if isuclid=false
if ~isuclid
    
    [lat1,lon1]=meshgrid(lat,lon);
    lat1=double(lat1);
    lon1=double(lon1);
    
    ngrid1=[0:dx:558];
    ngrid2=[0:dx:(366*(9/7))];
    [ngridx,ngridy] = meshgrid(ngrid1,ngrid2);
    
    
    lat=lat1; 
    lon=lon1;
    
    %Transpose obs data because interpolating requires 1st dimension be lon
    data=data'; 
    lat_ref=lat_ref';
    lon_ref=lon_ref';
    
    data_rs=reshape(data,size(data,1)*size(data,2),1);
    lat_rs=reshape(lat,size(lat,1)*size(lat,2),1);
    lon_rs=reshape(lon,size(lon,1)*size(lon,2),1);
%     lat_ref_rs=reshape(lat_ref,size(lat_ref,1)*size(lat_ref,2),1);
%     lon_ref_rs=reshape(lon_ref,size(lon_ref,1)*size(lon_ref,2),1);
%     grid1x_rs=reshape(cgridx,size(cgridx,1)*size(cgridx,2),1);
%     grid1y_rs=reshape(cgridy,size(cgridy,1)*size(cgridy,2),1);
    

    newdat=griddata(lat_rs,lon_rs,data_rs,lat_ref,lon_ref);
%     newlat=griddata(grid1x_rs,grid1y_rs,lat_rs,ngridx,ngridy);
%     newlon=griddata(grid1x_rs,grid1y_rs,lon_rs,ngridx,ngridy);
    newx=ngridx;
    newy=ngridy;
    %re transpose so we have interpolated data back in (lat,lon)
    %dimensional format
    newdat=newdat';
    
    
%else for isuclid=true    
else
    current_grid=dx_old; %grid spacing of orig grid
    interp_grid=dx; %new grid spacing
    %change 4/27/12 (...)-current_grid for cross-point grid
    cgrid1=[0:current_grid:(size(data,1)*current_grid)-current_grid];
    cgrid2=[0:current_grid:(size(data,2)*current_grid)-current_grid];

    [cgridy,cgridx] = meshgrid(cgrid2,cgrid1);
    
    newgrid_max1=ceil(floor(size(data,1)*current_grid/interp_grid)*interp_grid/current_grid);
    newgrid_max2=floor(ceil(size(data,2)*current_grid/interp_grid)*interp_grid/current_grid);
    
    dx1=size(data,1)*current_grid-newgrid_max1*current_grid;
    dx2=size(data,2)*current_grid-newgrid_max2*current_grid;

    %Uniquely create new grid going fom 3/7km to 9/7km spacing starting
    %from SW corner, cutting out ~0.4 km of top of model data for easier
    %interpolation of lat/lon
    ngrid1=[0:interp_grid:558];
    ngrid2=[0:interp_grid:(366*interp_grid)];

    [ngridy,ngridx] = meshgrid(ngrid2,ngrid1);

    %interpolate data onto new grid
    data_rs=reshape(data,1,size(data,1)*size(data,2));
%     lat_rs=reshape(lat,1,size(lat,1)*size(lat,2));
%     lon_rs=reshape(lon,1,size(lon,1)*size(lon,2));
    grid1x_rs=reshape(cgridx,1,size(cgridx,1)*size(cgridx,2));
    grid1y_rs=reshape(cgridy,1,size(cgridy,1)*size(cgridy,2));

    newdat=griddata(grid1x_rs,grid1y_rs,data_rs,ngridx,ngridy);
    
%     newlat=griddata(grid1x_rs,grid1y_rs,lat_rs,ngridx,ngridy);
%     newlon=griddata(grid1x_rs,grid1y_rs,lon_rs,ngridx,ngridy);
    newx=ngridx;
    newy=ngridy;
    
    %Trim unncessary column and row strips on N and E sides of domain
    newdat=newdat(1:434,1:366);
%     newlat=newlat(1:434,1:366);
%     newlon=newlon(1:434,1:366);
    
end


end

