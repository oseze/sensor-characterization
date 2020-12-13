%% Clean up
clear; close all; clc;
startup;
%% List all folders
% Finder: Cmd+K for at mounte bigNAS
inputPath = '/Volumes/ntch-optofluidic/2017/2017_oseze/LABT_20170508_oseze_NTCH_D0148'; 
%inputPath ='Z:\2017\2017_oseze\LABT_20170508_oseze_NTCH_D0148';
AllFiles = dir(fullfile(inputPath,'*.mat'));

Temp = struct('name','','CoM',[],'FWHM',[],'LFWHM',[],'maxIntensity',[]);
Boxes = struct('name','','comMap',[],'fwhmMap',[],'intMap',[],'LFWHM',[]);
cropped = struct('CoM',[],'X',[],'Y',[],'fwhm',[],'int',[],'name','','LFWHM',[]);

%% Import data
threshold = 0.5;
nFiles = numel(AllFiles);
for iFile = 1:nFiles
    Temp(iFile).name = AllFiles(iFile).name;
    Boxes(iFile).name = AllFiles(iFile).name;
    Data(iFile) = load(fullfile(inputPath,AllFiles(iFile).name)); %#ok
    Temp(iFile).maxIntensity = squeeze(max(Data(iFile).Raw.I))';
    [Temp(iFile).CoM,Temp(iFile).FWHM] = simpleCoM(Data(iFile).Raw.w',Data(iFile).Raw.I,threshold);
end
%% Fit FWHM with Lorentzian 
nLines = 100;
nFrames = 100;
%t0 = tic;
for iFile=1:numel(AllFiles)
    LFWHM=zeros(nFrames,nLines);
    %hWaitbar = waitbar(0,sprintf('Analyzing file %d',iFile));
    for iFrame=1:nFrames
        fprintf('File %d frame %d\n',iFile,iFrame)
        for iLine=1:nLines
            %waitbar(((iFrame-1)*nLines+iLine)/(nFrames*nLines));
            x=Data(iFile).Raw.w';
            y=Data(iFile).Raw.I(:,iLine,iFrame);
            startIntensity = 0;
            [startAmplitude,peakIdx] = max(y);
            startAmplitude = 0.05*startAmplitude;
            startWidth = 0.05;
            startPosition = x(peakIdx);
            fitResult = fitLorentzian(x,y,startIntensity,startAmplitude,...
                startWidth,startPosition,'PlotFit',false);
            
            LFWHM(iFrame,iLine) = abs(fitResult.width);         
        end        
    end
    saveAppendix = sprintf('_file%d',iFile);
    msave(saveAppendix,'LFWHM');
    Temp(iFile).LFWHM=LFWHM;
end
%tf = toc(t0)
msave('','Temp')

%%
dx = 0.5;
dy = 0.02;
[nx,ny]=size(Data(1).Results.map);
nx = nx/100;

%% Assemble maps 
for iFile = 1:nFiles
    iFrame = 1;
    Boxes(iFile).comMap = Data(iFile).Results.map;
    Boxes(iFile).x = Data(iFile).Results.x;
    Boxes(iFile).y = Data(iFile).Results.y;
    Boxes(iFile).intMap = 1000*ones(100*nx,ny);
    Boxes(iFile).fwhmMap = ones(100*nx,ny);
    Boxes(iFile).LFWHM = ones(100*nx,ny);

    for ix = 1:nx
        direction = -(-1)^ix;
        
        for iy = 1:ny
            xRange = 100*(ix-1)+1 : 100*ix;
             %map for CoM
            if direction > 0
                Boxes(iFile).comMap(xRange,iy) = Temp(iFile).CoM(iFrame,:);
            else
                Boxes(iFile).comMap(xRange,ny-iy+1) = Temp(iFile).CoM(iFrame,:);
            end   
            
            %map for int
            if direction > 0
                Boxes(iFile).intMap(xRange,iy) = Temp(iFile).maxIntensity(iFrame,:);
            else
                Boxes(iFile).intMap(xRange,ny-iy+1) = Temp(iFile).maxIntensity(iFrame,:);
            end
        
            
            %map for fwhm
            if direction > 0
                Boxes(iFile).fwhmMap(xRange,iy) = Temp(iFile).FWHM(iFrame,:);
            else
                Boxes(iFile).fwhmMap(xRange,ny-iy+1) = Temp(iFile).FWHM(iFrame,:);
            end
 
            %map for LFWHM
            if direction > 0
                Boxes(iFile).LFWHM(xRange,iy) = Temp(iFile).LFWHM(iFrame,:);
            else
                Boxes(iFile).LFHWM(xRange,ny-iy+1) = Temp(iFile).LFWHM(iFrame,:);
            end
            
            iFrame = iFrame + 1;
        end
    end
end

%% Remove the reference measurements
%%Boxes(1:2)=[];
%%
mFiles=numel(Boxes);
M=zeros(1,mFiles);
S=zeros(1,mFiles);
SFF=zeros(1,mFiles);
Mf=zeros(1,mFiles);
Mi=zeros(1,mFiles);
ML=zeros(1,mFiles);
SL=zeros(1,mFiles);

close all;
%% Crop image to make sure the meausurement are taken in the same place on the sensor
for iFile = 1:mFiles
close all
figure(iFile); imagesc(Boxes(iFile).comMap);
hColorbar = colorbar;
set(gcf,'OuterPosition',[1 58 1280 720]);

[x,y]= getpts;
%close(gcf)
p=4; h=80; w=20;

x=round(x+p); yh=round(y+h);
y=round(y+4*p); xw=round(x+w);

if yh > 200
    yh=200;
end

if xw >50
    xw=50;
end

cropped(iFile).CoM=Boxes(iFile).comMap(y:yh,x:xw); 
cropped(iFile).X=Boxes(iFile).x(y:yh); 
cropped(iFile).Y=Boxes(iFile).y(x:xw); 
cropped(iFile).fwhm=Boxes(iFile).fwhmMap(y:yh,x:xw); 
cropped(iFile).LFWHM=Boxes(iFile).LFWHM(y:yh,x:xw); 
cropped(iFile).int=Boxes(iFile).intMap(y:yh,x:xw); 
cropped(iFile).name=Boxes(iFile).name; 
cropped(iFile).LFWHM=Boxes(iFile).LFWHM; 
rectangle('Position',[x y xw-x yh-y]);
M(iFile)=mean2(cropped(iFile).CoM);
S(iFile)=std2(cropped(iFile).CoM);
SFF(iFile)=std2(cropped(iFile).fwhm);
Mf(iFile)=mean2(cropped(iFile).fwhm);
Mi(iFile)=mean2(cropped(iFile).int);
ML(iFile)=mean2(cropped(iFile).LFWHM);
SL(iFile)=std2(cropped(iFile).LFWHM);
figure; imagesc(cropped(iFile).CoM)
getpts;

end

saveAppendix2 = sprintf('_cropped');   
msave(saveAppendix2,'M','S','SFF','Mf','Mi','ML','SL','cropped')
%msave('','M','S','croppedCoM','croppedx','croppedy','croppedfwhm','croppedint','fname');