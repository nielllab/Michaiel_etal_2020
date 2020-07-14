clear all; close all;
dbstop if error
load('DEINTERLACED_multAniTest_071420.mat')
set(groot,'defaultFigureVisible','on') %disable figure plotting
deInter=1;


if deInter
frRate=60;
else
frRate=30; 
end

mouse_xy=[];cricket_xy=[];EllipseParamsR=[];EllipseParamsL=[]; radR=[]; az=[];

for i=1:length(Data)
    
    animalName(i,:)=Data(i).ani;
    expSession(i) = Data(i).sessionnum;
    expDate(i)=Data(i).date;
    clipNumber(i)=Data(i).clipnum;
    
     mouse_xy{i,1,1}=Data(i).mouseXY; 
%      mouseVel{i,:}= Data(i).mouseV; %in pix/frame
     mouseVel{i,:}=((Data(i).mouseV)/27)*frRate; %cm/sec, with our top camera 27 pix=1cm
     cricket_xy{i,1}=Data(i).cricketxy;
     cricketTheta{i,:}=Data(i).cricket_theta; 
%     %    cricketVel{i,:}= Data(i).cricketV % pix/frame
     cricketVel{i,:} = ((Data(i).cricketV)/27)*frRate; % cm/sec
     headTheta{i,:}= rad2deg(Data(i).theta);
     dTheta{i,:}=rad2deg(Data(i).dth);
     dist2cricket{i,:}= (Data(i).range)/27; %convert pixels to cm
    az = mod(Data(i).az,2*pi); az(az>pi) = az(az>pi)-2*pi;
    azT{i,:,:}= az; % as in radians
    azimuth{i,:,:}=rad2deg(azT{i,:,:});%az in deg
    slip{i,1} = Data(i).difTR;
    slip{i,2} = Data(i).difTL;
    slip{i,3} = Data(i).difRL;
    goodTheta(i,:) = Data(i).ThetaFract; %for each experiment, what is fraction of non-nan pts, from un-interpolated theta
    TSused{i,:}=Data(i).usedTS; %time stamps used for interpolation
    
    accelerometerChs{i,:} = Data(i).accShift;
    accelerometerChsRaw{i,:}= Data(i).rawAccShift;
    
    eyefitParamsR{i,1,1} = Data(i).EllipseParamsR;
    eyefitParamsL{i,1,1} = Data(i).EllipseParamsL;
    thetaR{i,:} = Data(i).Rtheta;
    thetaL{i,:} =Data(i).Ltheta;
    phiR{i,:} =Data(i).Rphi;
    phiL{i,:} =Data(i).Lphi;
   
    dthetaR{i,:} =(Data(i).dxRTheta);
    dthetaL{i,:} =(Data(i).dxLTheta);
    dphiR{i,:} =(Data(i).dxRPhi);
    dphiL{i,:} =(Data(i).dxLPhi);
     
%     XRcent{i,:} =Data(i).XRcent; %centroid px values of eyes, not really
%     used - use eye theta values instead
%     YRcent{i,:} =Data(i).YRcent;
%     XLcent{i,1} = Data(i).XLcent;
%     YLcent{i,:} =Data(i).YLcent;

%%% some checks on eye tracking and ellipse fit quality 
    goodR{:,i}=Data(i).goodReye; %1=all 8pts above likelihood .95
    Rngood{i,:}=Data(i).ngoodR; % num good DLC pts
    Rcc(i)=Data(i).RcalR; % corr coefficient for calibration values (calculated in EyeCameraCalc1, called in loadAllCsv)
    eyeSlopeR(i)=Data(i).RcalM; % slope for calibration values
    eyeScaleR(i)=Data(i).scaleR;
    rMissing(i)=sum(~isnan(Data(i).Rtheta))/(length(Data(i).Rtheta)); %proportion of non-nans, to filter out bad datasets    
%%% same as above but for L eye
    goodL{i,:}=Data(i).goodLeye;
    Lngood{i,:}=Data(i).ngoodL;
    Lcc(i)=Data(i).LcalR;
    eyeSlopeL(i)=Data(i).LcalM;
    eyeScaleL(i)=Data(i).scaleL;
    lMissing(i)=sum(~isnan(Data(i).Ltheta))/(length(Data(i).Ltheta));

end

%% select good datasets by setting thresholds on some QC measures
goodRight = Rcc>.3 & rMissing>.75; 
goodLeft=Lcc>.3 & lMissing>.75;

useTime = goodTheta>=.9 & goodRight' & goodLeft';
useData=find(useTime); 

allTilt=[];allRoll=[];allYaw=[]; dlcVerg=[];dlcDverg=[];allGyro1=[];allGyro2=[];allGyro3=[];
dlcDverg=[]; dlcDphi=[]; dlcPhi=[];dlcDhth=[];dlcHth=[];

for vid = 1:length(useData)
    roll = (accelerometerChs{useData(vid)}(:,1)); roll=roll-nanmean(roll);
    rollFilt = medfilt1(roll,2);
    
    tilt = (accelerometerChs{useData(vid)}(:,2)); tilt=tilt-nanmean(tilt);
    tiltFilt = medfilt1(tilt,8);
    
    yaw = (accelerometerChs{useData(vid)}(:,3)); yaw=yaw-nanmean(yaw);
    yawFilt = medfilt1(yaw,10);
    
    verg=(thetaR{useData(vid)}-thetaL{useData(vid)});
    dverg=(dthetaR{useData(vid)}-dthetaL{useData(vid)});
    
    gyro1=(accelerometerChs{useData(vid)}(:,4));
    gyro2=accelerometerChs{useData(vid)}(:,5);
    gyro3=(accelerometerChs{useData(vid)}(:,6));
    
    phi=(phiR{useData(vid)}-phiL{useData(vid)}); phi=phi-nanmean(phi);
    dphi=(dphiR{useData(vid)}-dphiL{useData(vid)}); dphi=dphi-nanmean(dphi);
    
    badE=dphi>15|dphi<-15 ; %noisy pts
    dphi(badE)=nan;
    
    hth=headTheta{useData(vid)}; hth=hth-nanmean(hth);
    dhth=dTheta{useData(vid)}; dhth=dhth-nanmean(dhth);
   
    bad=dhth>20 | dhth<-20; %noisy pts
    dhth(bad)=nan;
    
    allTilt=[allTilt tiltFilt(1:end-1)'];
    allRoll=[allRoll rollFilt(1:end-1)'];
    allYaw=[allYaw yawFilt(1:end-1)'];
    
    vergence=[dlcVerg verg(1:end-1)];
    dVergence=[dlcDverg dverg];
   
    allGyroRoll=[allGyro1 gyro1(1:end-1)']; 
    allGyroTilt=[allGyro2 gyro2(1:end-1)'];
    allGyroYaw=[allGyro3 gyro3(1:end-1)'];
    
    allAvgPhi=[dlcPhi phi(1:end-1)];
    d_allAvgPhi=[dlcDphi dphi];
    
    headThDLC=[dlcHth hth(1:end-1)];
    d_headThDLC=[dlcDhth dhth(1:end-1)'];
end
%% Identify Approach periods

for vid=1:length(useData)
    deltaR = diff(dist2cricket{useData(vid)})*frRate;
    vsmooth = conv(mouseVel{useData(vid)},ones(5,1)/5,'same');
    dRThresh=-10; %%% distance to cricket must be decreasing by 10 cm/sec
    vThresh=10; %%% speed of mouse must be at least cm/sec
    azThresh = pi/4;  %%% head angle relative to cricket must be pi/4 = 45 deg
    approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(azT{useData(vid)}(1:end-1))<azThresh;
    approach(1)=0; approach(end)=0; %%% boundary conditions
    
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop
    
    for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
        if (starts(j+1)-ends(j))<5 & (dist2cricket{useData(vid)}(starts(j+1))- dist2cricket{useData(vid)}(ends(j)))<3
            approach(ends(j) : starts(j+1))=1;
        end
    end
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop
    
    for j = 1:length(starts);   %%% remove short approaches (less than 10 frames)
        if ends(j)-starts(j)<10
            approach(starts(j):ends(j))=0;
        end
    end
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop
    approachEpochs{vid,:}=approach;
end

%% if you want to save out data from only good experiments

accelerometerChs=accelerometerChs(useData);
accelerometerChsRaw=accelerometerChsRaw(useData);
animalName=animalName(useData);
azimuth=azimuth(useData);
clipNumber=clipNumber(useData);
cricket_xy=cricket_xy(useData);
cricketVel=cricketVel(useData);
cricketTheta=cricketTheta(useData);
dTheta=dTheta(useData);
expDate=expDate(useData);
dist2cricket=dist2cricket(useData);
dphiL=dphiL(useData);
dthetaL=dthetaL(useData);
dphiR=dphiR(useData);
dthetaR=dthetaR(useData);
eyefitParamsL=eyefitParamsL(useData);
eyefitParamsR=eyefitParamsR(useData);
phiL=phiL(useData);
thetaL=thetaL(useData);
mouse_xy=mouse_xy(useData);
mouseVel=mouseVel(useData);
phiR=phiR(useData);
thetaR=thetaR(useData);
expSession=expSession(useData);
eyeScaleL=eyeScaleL(useData);
eyeScaleR=eyeScaleR(useData);
eyeSlopeL=eyeSlopeL(useData);
eyeSlopeR=eyeSlopeR(useData);
headTheta=headTheta(useData);



%%
pSname='/Users/angiemichaiel/Desktop/GitHubCode/Michaiel_etal_2020/';
afilename=('multAni_test_071420_COMPILED_b.mat');
save(fullfile(pSname, afilename))