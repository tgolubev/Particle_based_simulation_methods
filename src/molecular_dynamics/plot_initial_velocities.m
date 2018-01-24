%plot initial velocities

clear all

[File,Path]=uigetfile('*.txt','MultiSelect','off');


    str=sprintf('%s', [Path File]);   
     format shortG                                              %change formating so doesn't show 0's for e-11 values. 

     data = importdata(str);


     vx = data(:,1);
     vy = data(:,2); %the /10 is to convert A/m^2 to mA/cm^2
     vz = data(:,3);
   
 
        matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual') 
        
     scatter3(vx,vy,vz);   
     set(gca, 'FontSize', 24)
     xlabel('Position','interpreter','latex','FontSize',26.4);
     %ylabel({'Current Density ($mA/cm^2$)'},'interpreter','latex','FontSize',26.4)
     
     