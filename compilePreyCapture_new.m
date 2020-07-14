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
    
    animalName(i,:)=Data(i).ani;
    expSession(i) = Data(i).sessionnum;
    expDate(i)=Data(i).date;
    clipNumber(i)=Data(i).clipnum;
    
%      mouse_xy{i,1,1}=Data(i).mouseXY; -
%      mouseVel{i,:}= Data(i).mouseV; %in pix/frame
     mouseVel{i,:}=((Data(i).mouseV)/27)*frRate; %cm/sec, with our top camera 27 pix=1cm
     cricket_xy{i,1}=Data(i).cricketxy;
  %   cricketTheta{i,:}=Data(i).cricket_theta; -
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
    goodTheta(i,:) = Data(i).ThetaFract; %from un-interpolated theta
%     TSused{i,:}=Data(i).usedTS; %time stamps used for interpolation
%     thMissing(i)= sum(~isnan(Data(i).theta))/(length(Data(i).theta)); %proportion of non-nans
%     
%     accelData{i,:} = Data(i).accShift;
% %    accelDrift(i) =Data(i).accXcorrLag;
%     accelDataRaw{i,:}= Data(i).rawAccShift;
%     
%     EllipseParamsR{i,1,1} = Data(i).EllipseParamsR;
%     EllipseParamsL{i,1,1} = Data(i).EllipseParamsL;
%     thetaR{i,:} = Data(i).Rtheta;
%     thetaL{i,:} =Data(i).Ltheta;
%     phiR{i,:} =Data(i).Rphi;
%     phiL{i,:} =Data(i).Lphi;
%     rMissing(i)=sum(~isnan(Data(i).Rtheta))/(length(Data(i).Rtheta)); %proportion of non-nans
%     lMissing(i)=sum(~isnan(Data(i).Ltheta))/(length(Data(i).Ltheta));
%     
%     dthetaR{i,:} =(Data(i).dxRTheta);
%     dthetaL{i,:} =(Data(i).dxLTheta);
%     dphiR{i,:} =(Data(i).dxRPhi);
%     dphiL{i,:} =(Data(i).dxLPhi);
%     
%     dXRcent{i,:} =Data(i).dxR;
%     dXLcent{i,:} =Data(i).dxL;
%     %     dYRcent{i,:} =Data(i).dyR;
%     %     dYLcent{i,:} =Data(i).dyL;
%     XRcent{i,:} =Data(i).XRcent;
%     YRcent{i,:} =Data(i).YRcent;
%     XLcent{i,1} = Data(i).XLcent;
%     YLcent{i,:} =Data(i).YLcent;
%     RRad{i,:}=Data(i).RRad;
%     LRad{i,:}=Data(i).LRad;
% %     
% %     if deInter
% %     goodR{:,i}=Data(i).goodReye;
% %     else
% %     goodR{i,:}=Data(i).goodReye; %1=all 8pts above likelihood .95
% %     end
%     Rngood{i,:}=Data(i).ngoodR; % num good DLC pts
%     Rcc(i)=Data(i).RcalR;
%     Rslope(i)=Data(i).RcalM;
%     Rscale(i)=Data(i).scaleR;
%     
% %     goodL{i,:}=Data(i).goodLeye;
%     Lngood{i,:}=Data(i).ngoodL;
%     Lcc(i)=Data(i).LcalR;
%     Lslope(i)=Data(i).LcalM;
%     Lscale(i)=Data(i).scaleL;
%     
%     
%     %     longR{i,1}=EllipseParamsR{i,1}(:,3);longL{i,1}=EllipseParamsL{i,1}(:,3);
%     %     shortR{i,1}=EllipseParamsR{i,1}(:,4);shortL{i,1}=EllipseParamsL{i,1}(:,4);
%     %
%     %     radR{i,1}= (longR{i,1}+shortR{i,1})./length(longR{i,1});
%     %     radL{i,1}= (longL{i,1}+shortL{i,1})./length(longL{i,1});
%     %     pupilRvel{i,1}=diff(radR{i,1}); pupilLvel{i,1}=diff(radL{i,1})
%     
%     % tsData{i,1}=~isempty(Data(i).TopTs);
%     
end
%