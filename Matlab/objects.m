function [object_list,cluster_ct] = objects(x,y,cluster_ct,thresh,interpolated_grid)

%This function takes in a list of points (x and y) that have dBz greater than a
%pre-defined value (sorted) and groups points into objects and merges those
%that are within a spatial distance.
%Finds points within range of dBz point (going down x,y list), adds those
%as an object. Then functions loops through those points to find others
%within range to expand out the object.
%The list of x,y points will then have a corresponding number to what
%object it belongs.

object_list=zeros(length(x),1);
original_grid=3/7;
    
start=cluster_ct; %cluster_ct coming in from main script. Make sure > 0
if start==0
    start=1;
end

for i=1:length(x)
    if object_list(i)==0;
        cluster_ct=cluster_ct+1;
        object_list(i)=cluster_ct;
        %Calculate distances from i cell
        dx=(((x(i)-x(:)).^2)+((y(i)-y(:)).^2)).^0.5;
        neighbor_inds=find((dx<=thresh)&(object_list==0));
        
    if ~isempty(neighbor_inds) %if Not Empty
        %first set all inds to current cluster_ct
        object_list(neighbor_inds)=cluster_ct;
        %now loop through list of cells in current cluster_ct
        k=1; %prime loop
        while k < length(neighbor_inds)
            clear dx
            dx=(((x(neighbor_inds(k))-x(:)).^2)+((y(neighbor_inds(k))-y(:)).^2)).^0.5;
            secondary_inds=find((dx<=thresh)&(object_list==0));
            object_list(secondary_inds)=cluster_ct;
            neighbor_inds(length(neighbor_inds)+1:length(neighbor_inds)+length(secondary_inds))=secondary_inds;
            clear secondary_inds
            k=k+1;
        end        
    end
            %require object size be resolvable by grid scale
            cluster_list_inds=find(object_list==cluster_ct);
            
            if length(cluster_list_inds)*(interpolated_grid^2) < (original_grid*8)^2
                object_list(cluster_list_inds)=NaN;%set cluster_list_inds back to NaN!
                cluster_ct=cluster_ct-1; %set cluster_ct back to restart that integer value
            end
    clear neighbor_inds  
    end
end    
    %square fit min size of resolvable object by model
%     original_grid=3/7;
%     minsize=(original_grid*6)^2;
%     
%     for h=start:cluster_ct
%              %require object size be resolvable by grid scale
%             check_inds=find(object_list==h);
%             checksize=length(check_inds)*((interpolated_grid)^2);
%              if checksize < minsize
%                  object_list(check_inds)=NaN;
%               
%                  replace_inds=find(object_list>h);
%                  object_list(replace_inds)=object_list(replace_inds)-1;
%                  
%              end
%        clear check_inds
%        clear replace_inds
%     end
   cluster_ct=max(max(object_list));
   
   if isnan(cluster_ct)==1; %is equal to NaN
       cluster_ct=0; %all objects too small in list and set to NaNs
   end
   
end
