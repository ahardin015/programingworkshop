function [new] = new_object(rows, cols, refs, interpolated_grid, time, ob, t)


    new(1,1)=ob;
    new(1,2)=0; %CI flag (becomes 1 when matched 3 times)
    new(1,3)=0; %Matche count in time (tally)
    new(1,4)=str2double(time); %First time
    new(1,5)=str2double(time); %Latest time
    %Area average object center
    new(1,6)=sum(cols)/length(cols); %cols_center(ob)
    new(1,7)=sum(rows)/length(rows); %rows_center(ob)
    new(1,8)=(length(rows))*(interpolated_grid^2); %object_area(ob)
    new(1,9)=t;
end
