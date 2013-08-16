%%% Matlab Tutorial Script %%%
%%% you can either type the name of this script at the Matlab command line
%%% or click "run" above

%Comments are made using the "%" symbol 
%you can comment multiple lines by hightling and hit "Ctrl + r"
%to uncomment, highlight and hit "Ctrl + t" (Note: different for Macs)

clear all %this clears any variables being stored and is a good idea to have at the beginning of a script

%you can use Matlab as a calculator
%Matlab math syntax examples
5*5
5^2 %these are the same
125/5
5+5; %This answer will not be sent to the command line but rather stored as "ans" variable

%Mat lab's speciality is with matrix calculations
a=[1,2,3,4]; %this creates a column of cells, note ";" surpresses command line output
b=[4;3;2;1]; %this creates a row of cells

c=a*b; %1x4 array multiplied by 4x1 array gives 1x1 answer

%pre allocating will help with computation time when dealing with large
%matrices
d=zeros(100,100,10);

%how to write to command line
disp('milk was a bad choice')


x=rand(500,1); %generates 500 random values from 0-1

%plot it
figure %always start with this call
plot(x) %you can set plot multiple lines and with multiple dimensions
xlabel('This is your x-label') %x-axis label, make sure to use '' to create a string
ylabel('This is your y-label')
title('This is your title for the random number plot')
%to save the image
%print -dpdf random.pdf %this would save as a .pdf image 

%produce a histogram
figure
hist(x)

%basic statistics (click see Examples etc. at top of Command Window for
%examples, videos etc.)
max_value=max(x); %just max(x) will return the maximum value of every column
mean_value=mean(x);
 
%note: apply a function multiple times to matrices greater than 1-d
max_2d=max(max(d)); %finds max each column then finds the max out of those values
%%%


%producing map images (not the greatest)
figure
axesm lambert;
ax=usamap('texas');
latlim=getm(ax,'maplatlimit');
lonlim=getm(ax,'maplonlimit');
states=shaperead('usastatehi','UseGeoCoords',true,'boundingbox',[lonlim',latlim']);
geoshow(ax, states,'facecolor','none')




    
    
