clear all; close all;
dbstop if error
load('multAnimals_DEINTERLACED_061120_halfShift_a.mat')
set(groot,'defaultFigureVisible','on') %disable figure plotting
deInter=1;

if deInter
frRate=60;
else
frRate=30; 
end

mouse_xy=[];cricket_xy=[];EllipseParamsR=[];EllipseParamsL=[]; radR=[]; az=[];

for i=1:length(Data)
    
    animal(i,:)=Data(i).ani;
    sessionN(i) = Data(i).sessionnum;
    expdate(i)=Data(i).date;
    clip(i)=Data(i).clipnum;
    
    mouse_xyRaw{i,1,1}=Data(i).mouse_xyRaw;
    mouseVRaw{i,:}= Data(i).mouseVRaw; %in pix/frame
    mouseVRaw{i,:}=((Data(i).mouseVRaw)/27)*frRate; %cm/sec
    cricket_xyRaw{i,1}=Data(i).cricketxyRaw;
    cricketVRaw{i,:}= Data(i).cricketVRaw % pix/frame
    cricketVRaw{i,:} = ((Data(i).cricketVRaw)/27)*frRate; %now cm/sec
    thetaRaw{i,:}= Data(i).thetaRaw;
    dThetaRaw{i,:}= diff(Data(i).thetaRaw);
    rangeRaw{i,:}= (Data(i).rangeRaw)/27; %convert pixels to cm
    azRaw = mod(Data(i).azRaw,2*pi); azRaw(azRaw>pi) = azRaw(azRaw>pi)-2*pi;
    azTRaw{i,:,:}= azRaw;
    azDegRaw{i,:,:}=rad2deg(azTRaw{i,:,:});%az in deg
    %     crickHRaw{i,1,1} = Data(i).cricketHead;
    %     crickH{i,1,1} = Data(i).crickH;
    
    mouse_xy{i,1,1}=Data(i).mouse_xy;
    % mouseV{i,:}= Data(i).mouseV; %in pix/frame
    mouseV{i,:}=((Data(i).mouseV)/27)*frRate; %cm/sec
    cricket_xy{i,1}=Data(i).cricketxy;
    cricketAz{i,:}=Data(i).cricketTheta;
    %    cricketV{i,:}= Data(i).cricketV % pix/frame
    cricketV{i,:} = ((Data(i).cricketV)/27)*frRate; %now cm/sec
    theta{i,:}= rad2deg(Data(i).theta);
    dTheta{i,:}=rad2deg(Data(i).dth);
%     dThetaFilt{i,:}=medfilt1(dTheta{i,:},3);
    range{i,:}= (Data(i).range)/27; %convert pixels to cm
    az = mod(Data(i).az,2*pi); az(az>pi) = az(az>pi)-2*pi;
    azT{i,:,:}= az;
    azDeg{i,:,:}=rad2deg(azT{i,:,:});%az in deg
    slip{i,1} = Data(i).difTR;
    slip{i,2} = Data(i).difTL;
    slip{i,3} = Data(i).difRL;
    goodTheta(i,:) = Data(i).ThetaFract; %from un-interpolated theta
    TSused{i,:}=Data(i).usedTS;
    thMissing(i)= sum(~isnan(Data(i).theta))/(length(Data(i).theta)); %proportion of non-nans
    
    accelData{i,:} = Data(i).accShift;
  %  accelCorr(i) =Data(i).accXcorrMax;
%     accelDrift(i) =Data(i).accXcorrLag;
    accelDataRaw{i,:}= Data(i).rawAccShift;
    
    EllipseParamsR{i,1,1} = Data(i).EllipseParamsR;
    EllipseParamsL{i,1,1} = Data(i).EllipseParamsL;
    thetaR{i,:} = Data(i).Rtheta;
    thetaL{i,:} =Data(i).Ltheta;
    phiR{i,:} =Data(i).Rphi;
    phiL{i,:} =Data(i).Lphi;
    rMissing(i)=sum(~isnan(Data(i).Rtheta))/(length(Data(i).Rtheta)); %proportion of non-nans
    lMissing(i)=sum(~isnan(Data(i).Ltheta))/(length(Data(i).Ltheta));
    
    dthetaR{i,:} =(Data(i).dxRTheta);
    dthetaL{i,:} =(Data(i).dxLTheta);
    dphiR{i,:} =(Data(i).dxRPhi);
    dphiL{i,:} =(Data(i).dxLPhi);
    
    dXRcent{i,:} =Data(i).dxR;
    dXLcent{i,:} =Data(i).dxL;
    %     dYRcent{i,:} =Data(i).dyR;
    %     dYLcent{i,:} =Data(i).dyL;
    XRcent{i,:} =Data(i).XRcent;
    YRcent{i,:} =Data(i).YRcent;
    XLcent{i,1} = Data(i).XLcent;
    YLcent{i,:} =Data(i).YLcent;
    RRad{i,:}=Data(i).RRad;
    LRad{i,:}=Data(i).LRad;
%     
%     if deInter
%     goodR{:,i}=Data(i).goodReye;
%     else
%     goodR{i,:}=Data(i).goodReye; %1=all 8pts above likelihood .95
%     end
    Rngood{i,:}=Data(i).ngoodR; % num good DLC pts
    Rcc(i)=Data(i).RcalR;
    Rslope(i)=Data(i).RcalM;
    Rscale(i)=Data(i).scaleR;
    
%     goodL{i,:}=Data(i).goodLeye;
    Lngood{i,:}=Data(i).ngoodL;
    Lcc(i)=Data(i).LcalR;
    Lslope(i)=Data(i).LcalM;
    Lscale(i)=Data(i).scaleL;
    
    
    %     longR{i,1}=EllipseParamsR{i,1}(:,3);longL{i,1}=EllipseParamsL{i,1}(:,3);
    %     shortR{i,1}=EllipseParamsR{i,1}(:,4);shortL{i,1}=EllipseParamsL{i,1}(:,4);
    %
    %     radR{i,1}= (longR{i,1}+shortR{i,1})./length(longR{i,1});
    %     radL{i,1}= (longL{i,1}+shortL{i,1})./length(longL{i,1});
    %     pupilRvel{i,1}=diff(radR{i,1}); pupilLvel{i,1}=diff(radL{i,1})
    
    % tsData{i,1}=~isempty(Data(i).TopTs);
    
end
%
% tsData= cell2mat(tsData)
% delayFull=cell2mat(slip);
goodR= Rcc>.3 & rMissing>.75;
goodL=Lcc>.3 & lMissing>.75;


useTime = goodTheta>=.9 & goodR' & goodL' & thMissing'>.90; %tsData==1; %|(useL & useR)
useFilt=find(useTime); %useFilt=useFilt(1:4,6:end);

%%
allTilt=[];allRoll=[];allYaw=[]; dlcVerg=[];dlcDverg=[];allGyro1=[];allGyro2=[];allGyro3=[];
dlcDverg=[]; dlcDphi=[]; dlcPhi=[];dlcDhth=[];dlcHth=[];
% figure
for vid = 1:length(useFilt)
    figure('units','normalized','outerposition',[0 0 1 1])
    roll = (accelData{useFilt(vid)}(:,1)); roll=roll-nanmean(roll);
    rollFilt = medfilt1(roll,2);
    tilt = (accelData{useFilt(vid)}(:,2)); tilt=tilt-nanmean(tilt);
    tiltFilt = medfilt1(tilt,8);
    yaw = (accelData{useFilt(vid)}(:,3)); yaw=yaw-nanmean(yaw);
    yawFilt = medfilt1(yaw,10);
    verg=(thetaR{useFilt(vid)}-thetaL{useFilt(vid)});
    dverg=(dthetaR{useFilt(vid)}-dthetaL{useFilt(vid)});
    gyro1=(accelData{useFilt(vid)}(:,4));
    gyro2=accelData{useFilt(vid)}(:,5);
    gyro3=(accelData{useFilt(vid)}(:,6));
    phi=(phiR{useFilt(vid)}-phiL{useFilt(vid)}); phi=phi-nanmean(phi);
    dphi=(dphiR{useFilt(vid)}-dphiL{useFilt(vid)}); dphi=dphi-nanmean(dphi);
    badE=dphi>15| dphi<-15 ;
    dphi(badE)=nan;
    hth=theta{useFilt(vid)}; hth=hth-nanmean(hth);
    dhth=dTheta{useFilt(vid)}; dhth=dhth-nanmean(dhth);
    bad=dhth>20 | dhth<-20;
    dhth(bad)=nan;

    
    allTilt=[allTilt tiltFilt(1:end-1)'];
    allRoll=[allRoll rollFilt(1:end-1)'];
    allYaw=[allYaw yawFilt(1:end-1)'];
    dlcVerg=[dlcVerg verg(1:end-1)];
    dlcDverg=[dlcDverg dverg];
    allGyro1=[allGyro1 gyro1(1:end-1)'];
    allGyro2=[allGyro2 gyro2(1:end-1)'];
    allGyro3=[allGyro3 gyro3(1:end-1)'];
    dlcPhi=[dlcPhi phi(1:end-1)];
    dlcDphi=[dlcDphi dphi];
    dlcHth=[dlcHth hth(1:end-1)];
    dlcDhth=[dlcDhth dhth(1:end-1)'];
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end



%% identify approach!!!

for vid=1:length(useFilt)
    deltaR = diff(range{useFilt(vid)})*frRate;
    vsmooth = conv(mouseV{useFilt(vid)},ones(5,1)/5,'same');
    dRThresh=-10; %%%cm/sec
    vThresh=10;
    azThresh = pi/4;  %%% pi/4 = 45 deg
    approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(azT{useFilt(vid)}(1:end-1))<azThresh;
    approach(1)=0; approach(end)=0; %%% boundary conditions
    
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop
    
    for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
        if (starts(j+1)-ends(j))<5 & (range{useFilt(vid)}(starts(j+1))- range{useFilt(vid)}(ends(j)))<3
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
    appEpoch{vid,:}=approach;
end

%%
pSname='T:\PreyCaptureAnalysis\Data\';
afilename='multAni_test_071320.mat');
save(fullfile(pSname, afilename))