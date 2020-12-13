% Clean up
%clear; 
close all; clc;
startup;

% Load data
mload

%% Lav liste
liste=repmat([1 2 ],1,4);
%% Plot data - CoM data
mapType = 'CoM';
medium={'Air','Water'};

for iFile=1:numel(cropped)
    figure(iFile)
    imagesc(cropped(iFile).X,cropped(iFile).Y,cropped(iFile).CoM);
    
    cropped(iFile).name=cropped(iFile).name(1:end-4); %#ok
    
    if strcmp('Sucrose', medium{liste(iFile)}) && strcmp('Water', medium{liste(iFile-1)})
        diff(iFile)=M(iFile)-M(iFile-1); %#ok
        h2=title(sprintf('%s: %s: %s: %.02f nm: %.02f nm', mapType, cropped(iFile).name...
            , medium{liste(iFile)}, M(iFile),diff(iFile)),'fontsize',14);
        
    else
        h2=title(sprintf('%s: %s: %s: %.02f nm ', mapType, cropped(iFile).name...
            , medium{liste(iFile)}, M(iFile)),'fontsize',14);
    end
    
    set(h2,'Interpreter','none')
    colormap(morgenstemning);
    set(gca,'YDir','normal')
    xlabel('x [mm]');
    ylabel('y [mm]');
    hColorbar = colorbar;
    
    appendString=sprintf('%s_%s_%s', mapType, cropped(iFile).name, medium{liste(iFile)});
    mexport(appendString)
end

%% Plot fwhm data 
mapType = 'FWHM';
medium={'Air','Water'};

for iFile=1:numel(cropped)
    figure(iFile)
    imagesc(cropped(iFile).X,cropped(iFile).Y,cropped(iFile).fwhm);
    
    cropped(iFile).name=cropped(iFile).name(1:end-4); %#ok
    caxis([0 4]);
    
        h2=title(sprintf('%s: %s: %s: %.02f nm ', mapType, cropped(iFile).name...
            , medium{liste(iFile)}, Mf(iFile)),'fontsize',14);
  
    
    set(h2,'Interpreter','none')
    colormap(morgenstemning);
    set(gca,'YDir','normal')
    xlabel('x [mm]');
    ylabel('y [mm]');
    hColorbar = colorbar;
    
    appendString=sprintf('%s_%s_%s', mapType, cropped(iFile).name, medium{liste(iFile)});
    mexport(appendString)
end

%% Plot int data 
mapType = 'Intensity';
medium={'Air','Water'};

for iFile=1:numel(cropped)
    figure(iFile)
    imagesc(cropped(iFile).X,cropped(iFile).Y,cropped(iFile).int);
    
    cropped(iFile).name=cropped(iFile).name(1:end-4); %#ok
    
        h2=title(sprintf('%s: %s: %s: %.02f nm ', mapType, cropped(iFile).name...
            , medium{liste(iFile)}, Mi(iFile)),'fontsize',14);

    
    set(h2,'Interpreter','none')
    colormap(morgenstemning);
    set(gca,'YDir','normal')
    xlabel('x [mm]');
    ylabel('y [mm]');
    hColorbar = colorbar;
    
    appendString=sprintf('%s_%s_%s', mapType, cropped(iFile).name, medium{liste(iFile)});
    mexport(appendString)
end
%%
vec=[2,4,6,8];
M_w=M(vec)'; Mf_w=Mf(vec)'; Mi_w=Mi(vec)'; S_w=S(vec)';
M_L=ML(vec)'; S_L=SL(vec)'; SFF=SFF(vec)';
sens=M_w-M(vec-1)';

%'M','S','Mf','Mi','ML','SL'