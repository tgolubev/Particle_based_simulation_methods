%plot initial positions

clear all

[File,Path]=uigetfile('*.txt','MultiSelect','off');


    str=sprintf('%s', [Path File]);   
     format shortG                                              %change formating so doesn't show 0's for e-11 values. 

     data = importdata(str);


     x = data(:,1);
     y = data(:,2); %the /10 is to convert A/m^2 to mA/cm^2
     z = data(:,3);
   
 
        matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual') 
        
     scatter3(x,y,z);   
     set(gca, 'FontSize', 24)
     xlabel('Position','interpreter','latex','FontSize',26.4);
     %ylabel({'Current Density ($mA/cm^2$)'},'interpreter','latex','FontSize',26.4)
     
     
     
 