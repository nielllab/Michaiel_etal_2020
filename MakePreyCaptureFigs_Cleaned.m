close all;
clear all;
%load('Michaiel_et_al.,2020_fullDataset.mat')
load('multAni_test_071420_COMPILED_b.mat')
load('compileAllAnimals_CONTROL_090219.mat','mouseSp','appEpoch','thetaHead');
frRate=60;
accelData=1;
%% FIGURE 1: a freely moving eye and head tracking system

%%% Figure 1C: some variables computed for a single example trial in the
%%%dataset
% for vid=105
%     use=approachEpochs{vid};
%     range=2*1148:2*1655
%     tR=thetaR{vid}(:,range)-nanmedian(thetaR{vid}(:,range));
%     tL=thetaL{vid}(:,range)-nanmedian(thetaL{vid}(:,range));
%     tR=interpNan(tR,10,'linear'); tL=interpNan(tL,10,'linear');
%     mnEye=.5*(tR+tL); mnEye=interpNan(mnEye,10,'linear');
%     headTh=headTheta{vid}(:,range)-nanmedian(headTheta{vid}(:,range));
%     headTh=headTheta{vid}(:,range)-nanmedian(headTheta{vid}(:,range));
%     eyeVelocity=interpNan(diff(mnEye), 10,'linear');
%     
%     figure
%     subplot(7,1,1)
%     plot(tR); hold on
%     plot(tL);
%     plot(mnEye-nanmedian(mnEye));xlim([1 length(range)])
%     ylim([-32 32]);
%     title('horiz. eye position')
%     subplot(7,1,2)
%     plot(eyeVelocity);
%     xlim([1 length(range)]);
%     title('eye velocity')
%     subplot(7,1,3)
%     plot(dist2cricket{vid}(:,range));xlim([1 length(range)]); hold on
%     ylim([0 40])
%     plot([296 296],[-60 60],'g')
%     title('distance to cricket')
%     subplot(7,1,4)
%     plot(rad2deg(azimuth{vid}(:,range)));xlim([1 length(range)])
%     title('azimuth')
%     subplot(7,1,5);
%     plot(mouseVel{vid}(:,range))
%     title('speed')
%     xlim([1 length(range)])
%     subplot(7,1,6)
%     plot(headTh);
%     xlim([1 length(range)])
%     title('head yaw')
%     subplot(7,1,7)
%     plot(accelerometerChs{vid}(range,2)-nanmedian(accelerometerChs{vid}(range,2)));
%     xlim([1 length(range)])
%     title('head pitch')
% end

%%% Figure 1D: mean velocity control & w/cameras

for vid=1:length(thetaHead)
    prop(vid)=sum(~isnan(thetaHead{vid}))./length(thetaHead{vid});
end
useData=find(prop>.85); %85% of points needs to be present for exp to be used

for vid=1:length(appEpoch)
    appTime=appEpoch{vid};
    runningSmooth =medfilt1(mouseSp{vid},(frRate/4)); %smooth running speed in 250 ms window
    Cntrl_speed(vid,1) =sum(runningSmooth(:,appTime==0)>5)./length(runningSmooth(:,appTime==0));
    Cntrl_speed(vid,3) =nanmean(runningSmooth(:,appTime==0));
    
    if sum(appTime==1)>frRate
        Cntrl_speed(vid,2) =sum(runningSmooth(:,appTime==1)>5)./length(runningSmooth(:,appTime==1));
        Cntrl_speed(vid,4) =nanmean(runningSmooth(:,appTime==1));
        
    else
        Cntrl_speed(vid,2)=nan;
        Cntrl_speed(vid,4)=nan;
        
    end
end
clear mouseSp appEpoch useData appTime thetaHead vid runningSmooth %clear control vars

for vid=1:length(approachEpochs)
    appTime=approachEpochs{vid};
    w=gausswin(frRate/2); %smooth running speed in 250 ms window
    runningSmooth =medfilt1(mouseVel{vid}(:,1:end-1),frRate/4);
    
    Cam_speed(vid,1) = sum(runningSmooth(appTime==0)>5)./(length(runningSmooth(appTime==0)));
    
    
    Cam_speed(vid,3) = nanmean(runningSmooth(appTime==0));
    if sum(appTime==1)>frRate
        Cam_speed(vid,2) =sum(runningSmooth(appTime==1)>5)./(length(runningSmooth(appTime==1)));
        Cam_speed(vid,4) = nanmean(runningSmooth(appTime==1));
    else
        Cam_speed(vid,2)=nan;
        Cam_speed(vid,4)=nan;
    end
end

Cntrl_speed(Cntrl_speed(:,4)>60,4)=NaN; % noisy control data, speeds at a few timepoints go above 150 cm/sec - bad DLC tracking @ those times
Cntrl_err= nanstd(Cntrl_speed)/sqrt(length(approachEpochs));
Cam_err= nanstd(Cam_speed)/sqrt(length(approachEpochs));

figure;plot(2,Cam_speed(:,3),'ko'); hold on
plot(2,nanmean(Cam_speed(:,3)),'b*','Markersize',15);
plot(2,nanmean(Cam_speed(:,3))+(Cam_err(:,3)),'r*','Markersize',10);
plot(2,nanmean(Cam_speed(:,3))-(Cam_err(:,3)),'r*','Markersize',10);
plot(3,Cam_speed(:,4),'go')
plot(3,nanmean(Cam_speed(:,4)),'b*','Markersize',15)
plot(3,nanmean(Cam_speed(:,4))+(Cam_err(:,4)),'r*','Markersize',10);
plot(3,nanmean(Cam_speed(:,4))-(Cam_err(:,4)),'r*','Markersize',10);

plot(5,Cntrl_speed(:,3),'ko')
plot(5,nanmean(Cntrl_speed(:,3)),'b*','Markersize',15)
plot(5,nanmean(Cntrl_speed(:,3))+(Cntrl_err(:,3)),'r*','Markersize',10);
plot(5,nanmean(Cntrl_speed(:,3))-(Cntrl_err(:,3)),'r*','Markersize',10);

plot(6,Cntrl_speed(:,4),'go')
plot(6,nanmean(Cntrl_speed(:,4)),'b*','Markersize',15)
plot(6,nanmean(Cntrl_speed(:,4))+(Cntrl_err(:,4)),'r*','Markersize',10);
plot(6,nanmean(Cntrl_speed(:,4))-(Cntrl_err(:,4)),'r*','Markersize',10);
xlim([1 7]); axis square
legend('cam','cam app','control','control app');
title('Figure 1D: avg speed w/cameras and without for approach and non-approaches')

%%% Figure 1E: average number of captures per session 
%%%avg numbers of crickets for each animal, across sessions, from
%%%experiment spreadsheet recorded by experimenters
exp=[4.6	6.2	5.2	4.333333333	3.166666667	4.166666667	5.2	6.6	3.5	4.166666667];
cntrl=[4.933333333	6.366666667	5.366666667	5.666666667	5.833333333	4.5	5	6.6	5.166666667	5.333333333];
expM=mean(exp); errExp=std(exp)./(length(sqrt(exp)));
mnC=mean(cntrl); errC=std(cntrl)./(length(sqrt(cntrl)));
comp=[expM mnC];err=[errExp,errC];
figure
barweb(comp,err)
title('Figure 1E: Average num. of captures/session')
legend('camera','control')
%% FIGURE 2: Coordination of eyes during free movement

% Figure 2B: scatter plots of R & L eyes yaw vs pitch, approach and non-approach
%%
figure
for vid=1:length(approachEpochs)
    nframe = min(length(thetaR{vid}),length(thetaL{vid}));
    pR=phiR{vid}(1:nframe); pL=phiL{vid}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=thetaR{vid}(1:nframe); tL=thetaL{vid}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    clear use
    useN = find(approachEpochs{vid}==0); use = find(approachEpochs{vid}==1);
   
    subplot(1,2,1)
    plot(tR(useN(1:480:length(useN))),pR(useN(1:480:length(useN))),'.b'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    plot(tR(use(1:480:length(use))),pR(use(1:480:length(use))),'.g'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    title('R eye')
    xlabel('yaw (deg)'); ylabel('pitch(deg)');
    subplot(1,2,2)
    plot(tL(useN(1:480:length(useN))),pL(useN(1:480:length(useN))),'.b'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    plot(tL(use(1:480:length(use))),pL(use(1:480:length(use))),'.g'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    xlabel('yaw (deg)'); ylabel('pitch(deg)');
    title('L eye')
end

suptitle('Figure 2B: R & L eye yaw & pitch');

%%
% Figure 2C: overlaid trace of two eyes converging and diverging

range=295:1495; %20 seconds
figure
appTimes=find(approachEpochs{24}(:,range)==1);
rT=thetaR{24}(:,range);
lT=thetaL{24}(:,range);
plot(rT-nanmean(rT)); hold on;
plot(lT-nanmean(lT)); hold on;
xlim([1 length(range)])
legend('r theta','l theta');
title('Figure 2C: example trace of horiz. eye position')
%running trace corresponding to time period above
figure
w = gausswin(10);
y = filter(w,1,mouseVel{24}(:,range));
plot(y/8); hold on;
xlim([1 length(range)])
title('Figure 2C: example trace of running velocity')


gyroBias=nanmedian(allGyroYaw); % subtract off median of gyro, which is a small bias in measurement
allGyroYaw = allGyroYaw-gyroBias;

% Figure 2E: correlation of eye thetas & eye phis
clear corr lags corrAll corrAAll uselags uselagsA
corrAll=[]; corrAAll=[];
for vid=1:length(approachEpochs)
    nframe = min(length(dthetaR{vid}),length(dthetaL{vid}));
    dtR=dthetaR{vid}(1:nframe); dtL=dthetaL{vid}(1:nframe);
    nonapp=approachEpochs{vid}==0;
    use = (nonapp==1)';
    if sum(use)>3
        [corr lags]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
        uselags=(lags>=-frRate & lags<=frRate);
    else
    end
    use=approachEpochs{vid}==1;
    if sum(use)>3 & sum(~isnan(dtR(use)))>20
        [corrA lagsA]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
        uselagsA=(lagsA>=-frRate& lagsA<=frRate);
    else
    end
    
    if sum(uselags)==2*frRate+1 & sum(uselagsA)==2*frRate+1
        corrAll(vid,:)=corr(uselags); corrAAll(vid,:)=corrA(uselagsA);
    else
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
err= nanstd(corrAll)/(sqrt(length(approachEpochs))); errA = nanstd(corrAAll)/(sqrt(length(approachEpochs)));
shadedErrorBar(1:size(corrAll,2),nanmean(corrAll,1),err,'-b',1); hold on
shadedErrorBar(1:size(corrAAll,2),nanmean(corrAAll,1),errA,'-g',1); hold on
plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.7 .7]);
xlim([41 81]); 
ylim([-.2 .4])
axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');

legend(L,{'dTheta non-app','dTheta Approach',}); 
title('Figure 2E: mean between eye corr, change in horiz. position');


% Figure 2F: hist of difference in horizontal eye angle (aka vergence)
clear vergCounts vgDiff nBins vgHist
for c=0:1
    for vid=1:length(approachEpochs)
        nframe = min(length(thetaR{vid}),length(thetaL{vid}));
        nframe=min(nframe, length(approachEpochs{vid}));
        tR=thetaR{vid}(1:nframe); tL=thetaL{vid}(1:nframe);
        tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
        vg=tR-tL;
        useN= approachEpochs{vid}(1:nframe)==c;
        nBins=-90:5:90
        if sum(useN)>4
                vgHist = (hist(vg(useN),nBins))/sum(useN&~isnan(vg));
        else
            vgHist=NaN;
        end
        vgDiff(vid,:,c+1)=vgHist;
    end
end

vergCounts(:,1)=nansum(vgDiff(:,:,1),1)./(length(approachEpochs));
vergCounts(:,2)=nansum(vgDiff(:,:,2),1)./(length(approachEpochs));
vergErr=nanstd(vgDiff)./(sqrt(length(vgDiff))); vergErr=squeeze(vergErr);

figure;shadedErrorBar(nBins,vergCounts(:,1),vergErr(:,1),'b',1); hold on
shadedErrorBar(nBins,vergCounts(:,2),vergErr(:,2),'g',1); hold on
title('Figure 2F: vergence (R - L eye)'); axis square
xlabel('vergence (deg)'); ylabel('proportion of time');
ylim([0 .27])
trMn=squeeze(nanmean(vgDiff,2));


% Figure 2H: scatter of head pitch (aka tilt) and eye vergence
close all
figure
for vid=1:length(approachEpochs)
    clear tiltVid
    nframe = min(length(thetaR{vid}),length(thetaL{vid}));
    tR=thetaR{vid}(1:nframe); tL=thetaL{vid}(1:nframe);
    tR=nanmedian(tR)-tR; tL=nanmedian(tL)-tL;
    
    tiltVid=accelerometerChs{vid}(:,2);
    tiltVid = medfilt1(tiltVid,8);
    
    tiltVid=tiltVid-nanmedian(tiltVid);
    vg=tR-tL;
    clear use
    for c=0:1
        use = find(approachEpochs{vid}==c);
        if c==0
            plot(tiltVid(use(1:300:length(use))),vg(use(1:300:length(use))),'b.'); hold on;
            xlim([-60 60]);ylim([-60 60]);
            axis square
        else
            plot(tiltVid(use(1:300:length(use))),vg(use(1:300:length(use))),'g.');
        end
      
    end
end
title('Figure 2H: scatter of head pitch (aka tilt) and eye vergence')


clear tiltHist
for c=0:1
    for vid=1:length(approachEpochs)
        clear useN tilt tiltHist
        tilt=accelerometerChs{vid}(:,1);
        useN= approachEpochs{vid}==c;
        nframe=min(length(tilt), length(useN));
        useN=useN(1:nframe);tilt=tilt(1:nframe);
        tilt=tilt-nanmean(tilt);
        nBins=-90:5:90
        if sum(useN)>1
            tiltHist = (hist(tilt(useN),nBins))./sum((~isnan(tilt(useN))))% & find(useN)');
        else
            tiltHist=NaN;
        end
        tiltFull(vid,:,c+1)=tiltHist;
    end
end

tiltCounts(:,1)=nansum(tiltFull(:,:,1),1)./(length(approachEpochs));
tiltCounts(:,2)=nansum(tiltFull(:,:,2),1)./((length(approachEpochs))-sum(isnan(tiltFull(:,1,2))))
tiltErr=nanstd(tiltFull)./(sqrt(length(tiltFull))); tiltErr=squeeze(tiltErr);

figure; shadedErrorBar(nBins,tiltCounts(:,1),tiltErr(:,1),'b',1); hold on
shadedErrorBar(nBins,tiltCounts(:,2),tiltErr(:,2),'g',1); hold on
title('Figure 2H: pitch'); axis square
xlabel('pitch (deg)'); ylabel('proportion of time');
ylim([0 .27])

%% FIGURE 3: Coordination of eyes and head

% put together data from all experiments
d_mnEyeAll=[];allDLChead=[];mnEyeAll=[]; mouseSpAll=[];
for vid=1:length(approachEpochs)
    nframe=min(length(dthetaR{vid}),length(dthetaL{vid}));
    nframe=min(nframe, length(dTheta{vid}));
    mnEye =.5*(thetaR{vid}(1:nframe)+thetaL{vid}(1:nframe));
    mnEye=mnEye-nanmean(mnEye);
    mnEyeD =.5*(dthetaR{vid}(1:nframe)+dthetaL{vid}(1:nframe));
    mnEyeD=mnEyeD-nanmean(mnEyeD);
    dHead=dTheta{vid}(1:nframe); dHead=dHead-nanmean(dHead);
    allDLChead = [allDLChead dHead'];
    speed =mouseVel{vid}(1:nframe);
    mouseSpAll=[mouseSpAll speed];
    mnEyeAll=[mnEyeAll mnEye];
    d_mnEyeAll=[d_mnEyeAll mnEyeD];
end

clear corrYaw corrAllYaw err errA lagsYaw
for c=0:1
    for vid=1:length(approachEpochs)
        nframe = min(length(accelerometerChs{vid}(:,6)),length(dthetaR{vid}));
        nframe = min(nframe, length(dthetaL{vid}));
        nframe=min(nframe, length(approachEpochs{vid}));
        use=approachEpochs{vid}(1:nframe)==c;
        dth=dTheta{vid}(1:nframe);
        gyroyaw=accelerometerChs{vid}(1:nframe,6); dtR=dthetaR{vid}(1:nframe); dtL=dthetaL{vid}(1:nframe);
        mnEye=.5*(dtR+dtL); gyroyaw=gyroyaw-nanmean(gyroyaw); mnEye=mnEye-nanmean(mnEye);
        if sum(use)>5 & sum(~isnan(gyroyaw(use)))>10
            [corrYaw lagsYaw]= nanxcorr(gyroyaw(use),mnEye(use),frRate,'coeff');
            uselags=(lagsYaw>=-frRate& lags<=frRate);
        end
        if sum(uselags)==2*frRate+1
            corrAllYaw(vid,:,c+1)=corrYaw(uselags);
        else
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
err= nanstd(corrAllYaw(:,:,1))/(sqrt(length(corrAllYaw)));
errA= nanstd(corrAllYaw(:,:,2))/(sqrt(length(corrAllYaw)));
shadedErrorBar(1:size(corrAllYaw,2),nanmean(corrAllYaw(:,:,1),1),err,'-b',1); hold on
shadedErrorBar(1:size(corrAllYaw,2),nanmean(corrAllYaw(:,:,2),1),errA,'-g',1);

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]);
ylim([-.55 .3]);
% xlim([21 41]); %15 and 46 ==500 ms
% xlim([8.5 53.5])
xlim([frRate-((frRate/2)+1) frRate+((frRate/2)+1)])
xlabel('time'); ylabel('correlation coeff');
axis square
title('Figure 3B: change in head yaw and horiz. eye velocity correlation');

figure
for c=0:1
    clear use
    use = find(allAppT==c);
    plot(allGyroYaw(use(1:100:end)), d_mnEyeAll(use(1:100:end)),'.'); hold on
    axis([-15 15 -15 15]);
    axis square
    hold on
end
xlabel('head yaw (from gyro)'); ylabel('horiz. eye velocity');
title('Figure 3C: scatter of head yaw & mean horiz. eye velocity')

figure
clear nBins h use hp
allRunFilt=medfilt1(mouseSpAll,15); %filter running trace 
stationary = allRunFilt<1; %1 cm/sec threshold

nBins= -40:1:40;
hs=hist(d_mnEyeAll(stationary),nBins);
plot(nBins,hs/sum(stationary)); hold on; title('mn horiz. eye velocity, stationary'); axis square

hmv=hist(d_mnEyeAll(stationary==0),nBins);
plot(nBins,hmv/sum(stationary==0)); hold on; title('mn horiz. eye velocity'); axis square
legend('stationary','moving');
ylim([0 .9]); xlim([-10 10]);
title('Figure 3D: when mouse is stationary, so are eyes');

% Figure 3E: hist of head yaw app and non approach (not different)
figure
for c=0:1
    clear nBins h use
    use = find(allAppT==c);
    nBins=-30:2:30
    h=hist(allGyroYaw(use),nBins);
    headHist(:,c+1)=h;
    plot(nBins, h/length(allGyroYaw(use))); hold on; title('head yaw from gyro')
    title('Figure 3E: head angle is not diff between non-app and app')
end

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = find(allAppT==c);
    nBins= -60:2:60
    h=hist(mnEyeAll(use),nBins);
    binoc(:,c+1)=h;  
    
    subplot(1,2,1)
    plot(nBins,h/length(mnEyeAll(use))); hold on; title('3F: mn horiz. eye pos'); axis square
    plot([-20,-20],[.15,0],'k--'); %roughly binocular zone
    plot([20,20],[.15,0],'k--');
    
    nBins=-30:2:30
    hp=hist(d_mnEyeAll(use),nBins);
    subplot(1,2,2)
    plot(nBins, hp/length(d_mnEyeAll(use))); axis square; hold on; title('3G, panel 1: mn Eye velocity') %note this is between frames, multiply by framerate (60hz) to get velocity in deg/sec
    deltaBinoc(:,c+1)=hp;
end

% Figure 3G, panel 2:when head is still (+/- 15 deg/sec aka .25*60), what is eye
% velocity?

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = (allAppT==c);
    still=(allGyroYaw<.25 & allGyroYaw>-.25);
    nBins= -10:1:10;
    h=hist(d_mnEyeAll(use&still),nBins);
    stillHist(:,c+1)=h;
    plot(nBins,h/sum(use&still)); hold on; title('mn Eye Yaw'); axis square
    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    xlim([-10 10]);
    
    
    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    
    title('Figure 3G: when head angle is not changing, eyes are also still');
    
end

%% FIGURE 4: Non-compensatory eye movements

%%% clustering of eye movement

dGz= allGyroYaw(allAppT==1) + d_mnEyeAll(allAppT==1); % change in gaze = sum of head yaw & horiz eye velocity

%%% threshold on gaze velocity
dGzV = dGz*frRate;
gzThresh = 180;  %% 
clust1 = abs(dGzV)>gzThresh;


gzVall = (allGyroYaw + d_mnEyeAll)*frRate;
figure
[gzHist bins] = hist(gzVall,-600:5:600);
semilogy(bins,gzHist/nansum(gzHist));
xlabel('gaze velocity deg/sec'); hold on
sPts = bins>=gzThresh;
semilogy(bins(sPts),gzHist(sPts)/nansum(gzHist),'r');
sPts = bins<=-gzThresh;
semilogy(bins(sPts),gzHist(sPts)/nansum(gzHist),'r');
plot([gzThresh gzThresh],[0.0001 0.1],'r:')
plot(-[gzThresh gzThresh],[0.0001 0.1],'r:')
title('Figure 4A: gaze velocity distribution, with threshold')

figure
gscatter(allGyroYaw(allAppT==1)-gyroBias,d_mnEyeAll(allAppT==1),clust1); axis equal; hold on
title('Figure 4B: scatter plot of yaw and mn horiz. eye position, approach pts only')
xlim([-25 25]); ylim([-25 25]); hold on; plot([-25 25],[25 -25],'g')

clear mvmts
mvmts=[allGyroYaw(allAppT==1); d_mnEyeAll(allAppT==1)];
gm = fitgmdist(mvmts',3,'Replicates',10); %gaussian mixture to cluster
idx = cluster(gm,mvmts'); 
dGzFull= allGyroYaw + d_mnEyeAll;
fullData=[allGyroYaw; d_mnEyeAll];
idxAll=cluster(gm, fullData');
[sacc, clust]= max([nanmean(abs(dGzFull(idxAll==1))),nanmean(abs(dGzFull(idxAll==2))), nanmean(abs(dGzFull(idxAll==3)))])
%%
maxlag = frRate;
thbins = -60:5:60;
skip = 1; %%% only shows figures at this interval
nthresh  = 60;
headAll=[];
eyeAll=[]; g3All= []; azAll = []; spAll  = []; allAppT = []; dist2cricketAll = [];
ns = 0; saccHeadAll = []; saccEyeAll = []; saccallAppT = []; saccVidAll= []; saccAzAll = []; saccThAll = []; saccEyeRawAll=[]; timetoallAppT =[];
% figure

%%% non-approach example in paper is vid 98, approach example is vid 42
%%% Note: azimuth to cricket is often empty, these are usually times when
%%% cricket is occluded (under mouse for example)
allS = 0; clear dheadStable eyeStable dgzStable %%% store out stable periods
for vid = 1:length(approachEpochs)
  
    %%% get approaches
    app = approachEpochs{vid};
    nonapp=~approachEpochs{vid};
    %%% get left eye positions
    lth = thetaL{vid} - nanmedian(thetaL{vid});
    dlth = dthetaL{vid};
    nl(vid) = sum(~isnan(lth(app))); %%% # good eye approach points
    lthHist(:,1,vid) = hist(lth(app),thbins)/nl(vid);
    lthHist(:,2,vid) = hist(lth(~app),thbins)/sum(~isnan(lth(~app)));
    if nl(vid)<nthresh
        lthHist(:,:,vid) = NaN;
    end
    %%% get right eye positions
    rth = thetaR{vid} - nanmedian(thetaR{vid});
    drth = dthetaR{vid};
    nr(vid) =sum(~isnan(rth(app))); %%% # good eye approach points
    rthHist(:,1,vid) = hist(rth(app),thbins)/nr(vid);
    rthHist(:,2,vid) = hist(rth(~app),thbins)/sum(~isnan(rth(~app)));
    if nr(vid)<nthresh
        rthHist(:,:,vid) = NaN;
    end
    
    
    %%% get eye phi
    Rphi = phiR{vid}-nanmean(phiR{vid});
    d_phiR = dphiR{vid};
    nrp(vid) =sum(~isnan(Rphi(app))); %%% # good eye approach points
    phiRHist(:,1,vid) = hist(Rphi(app),thbins)/nrp(vid);
    phiRHist(:,2,vid) = hist(Rphi(~app),thbins)/sum(~isnan(Rphi(~app)));
    if nr(vid)<nthresh
        phiRHist(:,:,vid) = NaN;
    end
    
    % get left eye phi
    Lphi = phiL{vid}-nanmean(phiL{vid});
    dphiL = dPhiL{vid};
    nlp(vid) =sum(~isnan(Lphi(app))); %%% # good eye approach points
    phiLHist(:,1,vid) = hist(Lphi(app),thbins)/nrp(vid);
    phiLHist(:,2,vid) = hist(Lphi(~app),thbins)/sum(~isnan(Lphi(~app)));
    if nr(vid)<nthresh
        phiLHist(:,:,vid) = NaN;
    end
    
    
    %%% get head positions
    hth = headTheta{vid};
    dth = dTheta{vid};
    azdeg = azimuth{vid}*180/pi;
    crickD = dist2cricket{vid};
    
    %%% get accelerometers
  
    if exist('accelData','var')
        tilt = accelerometerChs{vid}(:,2);
        roll = accelerometerChs{vid}(:,1);
        acc_dth = accelerometerChs{vid}(:,6);
        acc_dth=acc_dth-gyroBias;
        acc_hth = nancumsum([hth(1) acc_dth(1:end-1)'],[],2);
    else
        display('no acc')
    end
    
    
    %%% get head positions
    %     hth = thetaHead{vid}; %head angle from DLC
    %     dth = d_headThDLC{vid};
    %     azdeg = az{vid}*180/pi;
    
    hth = acc_hth; %head angle from accelerometer
    dth = acc_dth;
    azdeg = azimuth{vid}*180/pi;
    
    %%% azimuth vs eye histograms
    az_hist(:,1,vid) = hist(-azdeg(app),thbins)/sum(~isnan(azdeg(app)));
    n= sum(~isnan(azdeg(app)) & ~isnan(lth(app)));
    azthL_hist(:,1,vid) = hist(-azdeg(app)-lth(app),thbins)/n;
    n= sum(~isnan(azdeg(app)) & ~isnan(rth(app)));
    azthR_hist(:,1,vid) = hist(-azdeg(app)-rth(app),thbins)/n;
    
    %%% alignment of eyes during approaches
    vergence = rth-lth;
    n= sum(~isnan(vergence(app)));
    vergeHist(:,1,vid) = hist(vergence(app),thbins)/n;
    vergeHist(:,2,vid) = hist(vergence(~app),thbins)/sum(~isnan(vergence(~app)));
    
    %%% mean eye theta is most important for stabilization
    mnEyeTh = 0.5*(rth+lth);
    n= sum(~isnan(azdeg(app)) & ~isnan(mnEyeTh(app)));
    azthRL_hist(:,1,vid) = hist(-azdeg(app)-mnEyeTh(app),thbins)/n;
    
    %%% gaze is the sum of head position + mean eye position
    
        gaze = hth + mnEyeTh;
   
    %%% find the longest approach, to use as example image
    appPts = find(app);
    newApp = [1 find(diff(appPts)>1)+1];
    endApp = [find(diff(appPts)>1)-1 length(appPts)];
    
    dur = endApp - newApp;
    [mx longest] = max(dur);
    mainApp = appPts(newApp(longest)  : endApp(longest)); %%% approach time only
    %%% add 2 secs on either side
    try
        appStart = max(appPts(newApp(longest))-(2*frRate),1); appOffset = appPts(newApp(longest))-appStart;
        appEnd = min(appPts(endApp(longest))+(2*frRate),length(app)); endOffset = appPts(endApp(longest))-appStart;
        appRange = appStart  : appEnd;
    catch
        appRange = appPts(newApp(longest)  : endApp(longest));
        appOffset =0;
    end
    
    
    %%% get rid of large jumps
    
    hthnonan = acc_hth;
    
    hthApp = acc_hth(appRange);
    hthApp = hthApp - nanmedian(acc_hth(mainApp));
    
    
    gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
    mnEyeApp =  0.5*(rth(appRange) +lth(appRange))';
    
    dhthnonan = acc_dth;
    dhthnonan(abs(diff(dth))>90)=NaN;
    
    dhthApp = acc_dth(appRange)-nanmedian(acc_dth(mainApp));
    dhthApp = mod(dhthApp + 180,360)-180;
    
    %%% calculate change in position at different lags, as measure of stability
    
    %%% draw figures
    
    if round(vid/skip)==vid/skip
      
        roll = accelerometerChs{vid}(:,1); roll=roll-nanmean(roll); roll=medfilt1(roll,5);
        tilt=accelerometerChs{vid}(:,2); tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
        g3=accelerometerChs{vid}(:,6); g3=g3-gyroBias;
        
        dGaze = g3(1:end-1) + 0.5*(drth+dlth)';
     
        %%% assign clusters based on previous GMM; clust for analysis is set
        %%% above
        
        newData=[g3(1:end-1)' ; 0.5*(drth+dlth)];

        
    %%% threshold on gaze velocity
        dGazeV = dGaze*frRate;
        idx2 = abs(dGazeV)>=gzThresh;
        clust = 1;
        
        %%% get saccade starts
        sacc = idx2==clust;
        saccStart = [0; diff(sacc)>0];
        saccs = find(saccStart);
        clear saccHead saccEye saccApp saccVid saccAz saccTh saccEyeRaw timetoApp
        
        %%% grab traces around saccades;
       
    
        buffer = 20; %if frame rate is 60hz, divide by half if fr rate is 30 
        azAll = [azAll; azdeg'];
        eyeAll = [eyeAll; mnEyeTh'];
        g3All = [g3All; g3];
        allAppT = [allAppT; [app 0]' ];
        spAll = [spAll; mouseVel{vid}'];
        dist2cricketAll = [dist2cricketAll; dist2cricket{vid}'];
        
        if length(mouseVel{vid}) ~= length(g3)
            keyboard
        end
        
        
        for i = 1:length(saccs);
            
            if saccs(i)>buffer & saccs(i)<length(drth)-buffer
                ns = ns+1;
                saccHead(:,ns) = cumsum(g3(saccs(i) + (-buffer:buffer)));
                saccEye(:,ns) = mnEyeTh(saccs(i) + (-buffer:buffer));
                saccAz(:,ns) = azdeg(saccs(i) + (-buffer:buffer));
                saccTh(:,ns) = hth(saccs(i) + (-buffer:buffer));
                saccEyeRaw(:,ns) = saccEye(:,ns); %%% copy that is not inverted
                if dGaze(saccs(i))<0
                    saccHead(:,ns) = - saccHead(:,ns);
                    saccEye(:,ns) = - saccEye(:,ns);
                end
                saccApp(ns) = app(saccs(i));
                saccVid(ns) = vid;
                distToApp = saccs(i)-find(app);
                [d nearApp] = min(abs(distToApp));
                if isempty(nearApp)
                    timetoApp(ns)=NaN;
                else
                    timetoApp(ns) = distToApp(nearApp);
                end
            end
        end
        if exist('saccEye','var')
            saccEyeAll = [saccEyeAll saccEye];
            saccHeadAll = [saccHeadAll saccHead];
            saccallAppT = [saccallAppT saccApp];
            saccVidAll = [saccVidAll saccVid];
            saccAzAll = [saccAzAll saccAz];
            saccThAll = [saccThAll saccTh];
            saccEyeRawAll = [saccEyeRawAll saccEyeRaw];
            timetoallAppT = [timetoallAppT timetoApp];
        end
        
        %%% calculate stabilization;
        saccEnds = find(diff(sacc)<0)+1;
        if length(saccEnds)>0 & length(saccs>0) &  (saccEnds(1) < saccs(1)) %%% sometimes trace starts with a saccade end, messes things up
            saccEnds = saccEnds(2:end);
        end
        
        for s = 1:length(saccs)-1;
            stable = saccEnds(s):saccs(s+1)-1;
            allS = allS+1;
            dheadStable{allS} = g3(stable);
            dgzStable{allS} = dGaze(stable);
            eyeStable{allS} = mnEyeTh(stable);
            %             if ~isnan(sum(dgzStable{allS})) & std(cumsum(dgzStable{allS}))>25
            %                 keyboard
            %             end
        end
        
        figure('Renderer', 'painters', 'Position', [100 -100 600 900])
        subplot(5,1,1)
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on;
        plot(rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        plot(find(idx2(appRange)==clust),0.5*(rth(appRange(idx2(appRange)==clust)) +lth(appRange(idx2(appRange)==clust)))','og')
        ylim([-30 30]); % xlim([0 max(length(appRange),1)]);
        title('horiz. eye position')
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('mean','right','left');
        title(sprintf('vid %d',vid));
        
        subplot(5,1,2)
        plot(diff(rth(appRange)),'r');hold on; plot(diff(lth(appRange)),'b');
        plot(diff(0.5*( rth(appRange) +lth(appRange)))','k','LineWidth',2); hold on;
        title('eye velocity')
        
        subplot(5,1,3)
        plot(azdeg(appRange)); hold on
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        title('az to cricket')
        xlim([0 max(length(appRange),1)]);
        
        gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
        subplot(5,1,4)
        hold on
        plot(hthApp,'Color',[0 0.75 0],'LineWidth',2);
        plot((gzApp'+30),'k','LineWidth',2);
        title('gaze with saccades, and head yaw (green)')
        sacc = find(idx2(appRange)==clust);
        sacc = sacc(sacc<length(gzApp)-1); %%% make sure we don't run off the end of the data
        for s = 1:length(sacc)
            plot(sacc(s):sacc(s)+1,gzApp(sacc(s):sacc(s)+1) + 30,'r','LineWidth',2)
        end
        
        if appOffset>0
            plot([appOffset appOffset],[-60 60],'g');end
        plot([endOffset endOffset],[-60 60],'r');
        if ~isempty(min(hthApp)) & sum(~isnan(hthApp))>10
            ylim([min(hthApp)-20 max(gzApp)+30]);
        else  ylim([-60 60]);end
        % xlim([0 max(length(appRange),1)]);
        xlabel('frames'); ylabel('deg');

        subplot(5,1,5)
        gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
        plot((gzApp'),'k','LineWidth',2); hold on;
        sacc = find(idx2(appRange)==clust);
        sacc = sacc(sacc<length(gzApp)-1); %%% make sure we don't run off the end of the data
        for s = 1:length(sacc)
            plot(sacc(s):sacc(s)+1,gzApp(sacc(s):sacc(s)+1),'r','LineWidth',2)
        end
        
        if appOffset>0
            plot([appOffset appOffset],[-60 60],'g');end
        plot([endOffset endOffset],[-60 60],'r');
        if ~isempty(min(hthApp))
            ylim([min(gzApp)-5 max(gzApp)+5]);
            %     ylim([-50 50])
        else  ylim([-50 50]);end
        title('gaze w/saccades')
    end
    drawnow
   
end

%% analyze eye/head position around saccades
%%% defining 3 clusters


%%% can either select all saccs for analysis, or just approaches.
%%% right now, all saccs
useS = saccallAppT>=0;
eyeData = saccEyeAll(:,useS);
headData = saccHeadAll(:,useS);


%% calculate stability and frequency of fixations

%%% remove stabilizations with NaNs (since these may cause to miss a saccade)
for i = 1:length(eyeStable)
    nanCount(i) = sum(isnan(eyeStable{i})) + sum(isnan(dheadStable{i}));
end
sum(nanCount==0)/sum(nanCount>=0);
good = find(nanCount==0);
for i = 1:length(good)
    dheadStableGood{i} = dheadStable{good(i)};
    eyeStableGood{i} = eyeStable{good(i)};
    dgzStableGood{i} = dgzStable{good(i)};
end

%%% calculate stats for each fixation (duration and stdev of positions)
clear stabDur headStd gazeStd eyeStd headMn gazeMn
for i = 1:length(dheadStableGood);
    stabDurGood(i) = length(dheadStableGood{i});
    if length(dheadStableGood{i})>1   %%% can't do stats on one point
        headStd(i) = std(cumsum(dheadStableGood{i}));
        gazeStd(i) = std(cumsum(dgzStableGood{i}));
        eyeStd(i) = std(eyeStableGood{i});
        headMn(i)= nanmean((dheadStableGood{i}));
        gazeMn(i) = nanmean((dgzStableGood{i}));
    else
        headStd(i) = NaN;
        gazeStd(i) = NaN;
        eyeStd(i) = NaN;
        headMn(i)= NaN;
        gazeMn(i)= NaN;
    end
end

%%% plot histogram of durations (with error bars)
clear hbins
hbins = 1:2:150;
h =hist(stabDurGood,hbins); 
figure
shadedErrorBar(hbins/30,h/sum(h),sqrt(h)/sum(h)); xlabel('secs'); ylabel('fraction');
hold on; plot([nanmedian(stabDurGood)/frRate nanmedian(stabDurGood)/frRate], [0 .14])
title(sprintf('Figure 4D: fixation duration median = %0.2f +/- %0.2f sec',nanmedian(stabDurGood)/frRate,(1/frRate)*std(stabDurGood)/sqrt(length(stabDurGood)))); xlim([0 2]);


%%% histogram of stability for head and gaze (with error bars)
clear hbins
figure
hbins = 0.25:0.5:20;
hold on
h = hist(headStd,hbins);  %%% add error bars here!
shadedErrorBar(hbins,h/sum(h),sqrt(h)/sum(h),'b');
h = hist(gazeStd,hbins);
shadedErrorBar(hbins,h/sum(h),sqrt(h)/sum(h),'k');
legend('head','gaze');
xlabel('RMS stabilization (deg)'); ylabel('fraction'); xlim([0 15])
title('Figure 4E: RMS distributions')

%%% calculate average stability of head and gaze
stability(1) = nanmedian(headStd);
stability(2) = nanmedian(gazeStd);
stabErr(1) = nanstd(headStd)/sqrt(sum(~isnan(headStd)));
stabErr(2) = nanstd(gazeStd)/sqrt(sum(~isnan(gazeStd)));
figure
bar(1:2, stability);
hold on
errorbar(1:2,stability, stabErr,'.');
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'head','gaze'}); 
ylabel('RMS stability (deg)'); title('Figure 4F: stability during fixations')
apps = find(saccallAppT==1);

s = saccAzAll(:,apps);
azOffset=nanmedian(s(:))

if nanmean(azAll)<=0
    azAll = azAll-azOffset;  %%% don't do it twice
end
eyeAzAll = azAll+eyeAll;

use = abs(azAll)<90 & abs(g3All)<1;
nanmean(abs(azAll(use)))
nanmean(abs(eyeAzAll(use)))

spAll(spAll>50)= NaN;

range = 1:10:length(azAll);
spSmooth = medfilt1(spAll,5);
azSmooth = medfilt1(azAll,5);
azData = azSmooth(range); spData = spSmooth(range); appData = allAppT(range);
distData = dist2cricketAll(range);


%% FIGURE 5: eye and head relative to cricket for saccades during approaches
apps = find(saccallAppT==1);

%%% there is ~6deg asymmetry in head position relative to cricket. this is
%%% probably based on the definition of theta from "model" head points.
s = saccAzAll(:,apps);
azOffset=nanmedian(s(:))
saccAzAllRaw = saccAzAll; %%% save raw value
saccAzAll = saccAzAll - azOffset;

%%% get eyes relative to cricket
eyeAz = saccAzAll + saccEyeRawAll;
%%
%%% head and eye relative to cricket (azimuth) on all saccades
figure
trange=3:34
hold on
data = nanmedian(abs(eyeAz(:,apps)),2);
err = nanstd(abs(eyeAz(:,apps)),[],2) ./ sqrt(sum(~isnan(eyeAz(:,apps)),2))
data = data(trange); err = err(trange);
t = (0:length(data)-1)/frRate;
shadedErrorBar(t,data,err,'k')
data = nanmedian(abs(saccAzAll(:,apps)),2);
err = nanstd(abs(saccAzAll(:,apps)),[],2) ./ sqrt(sum(~isnan(saccAzAll(:,apps)),2))
data = data(trange); err = err(trange);
t = (0:length(data)-1)/frRate;
shadedErrorBar(t,data,err,'b')
legend('head azimuth','gaze azimuth');
ylim([7.5 27.5]);
xlim([t(1)-1/(2*frRate) t(end)+1/(2*frRate)]);
plot([t(19) t(19)],[7 28])
axis square
xlabel('secs');
ylabel('azimuth to cricket (deg)')
title('Figure 5B: Azimuth of gaze and head yaw relative to cricket, at time of saccade')

%%
%%% plot for all 4 clusters

%%% where does eye end up after saccade?
preSaccT = 21; postSaccT=23; %%% timepoints to compare pre/post saccade
bins = -(2*frRate):10:(2*frRate);

%%% get histograms of eye&head for pre/post saccade
clear eyeAzHist
eyeAzHist(:,1) = hist(eyeAz(preSaccT,apps),bins);
eyeAzHist(:,2) = hist(eyeAz(postSaccT,apps),bins);

AzHist(:,1) = hist(saccAzAll(preSaccT,apps),bins);
AzHist(:,2) = hist(saccAzAll(postSaccT,apps),bins);

%%% plot pre (with error bars)
figure
shadedErrorBar(bins,eyeAzHist(:,1)/sum(eyeAzHist(:,1)), sqrt(eyeAzHist(:,1))/sum(eyeAzHist(:,1)),'b')
hold on
shadedErrorBar(bins,AzHist(:,1)/sum(AzHist(:,1)), sqrt(AzHist(:,1))/sum(AzHist(:,1)),'r')
xlabel('azimuth (deg)'); legend('gaze','head')
title('pre saccade')
xlim([-60 60])
title('Figure 5C: az. distribution before saccade')

%%% plot post (with error bars)
figure
shadedErrorBar(bins,eyeAzHist(:,2)/sum(eyeAzHist(:,2)), sqrt(eyeAzHist(:,2))/sum(eyeAzHist(:,2)),'b')
hold on
shadedErrorBar(bins,AzHist(:,2)/sum(AzHist(:,2)), sqrt(AzHist(:,2))/sum(AzHist(:,2)),'r')
xlabel('azimuth (deg)'); legend('gaze','head')
title('post saccade')
xlim([-60 60])
title('Figure 5C: az. distribution after saccade')


%%% get median offset of eye&head for pre/post saccade
eyeOffset(1) = nanmedian(abs(eyeAz(preSaccT,apps)));
eyeOffset(2) = nanmedian(abs(eyeAz(postSaccT,apps)))
eyeOffsetErr(1) = nanstd(abs(eyeAz(preSaccT,apps)))/sqrt(sum(~isnan(eyeAz(preSaccT,apps))));
eyeOffsetErr(2) = nanstd(abs(eyeAz(postSaccT,apps)))/sqrt(sum(~isnan(eyeAz(postSaccT,apps))))

headOffset(1) = nanmedian(abs(saccAzAll(preSaccT,apps)));
headOffset(2) = nanmedian(abs(saccAzAll(postSaccT,apps)));
headOffsetErr(1) = nanstd(abs(saccAzAll(preSaccT,apps)))/sqrt(sum(~isnan(saccAzAll(preSaccT,apps))));
headOffsetErr(2) = nanstd(abs(saccAzAll(postSaccT,apps)))/sqrt(sum(~isnan(saccAzAll(postSaccT,apps))))

eyez=nanmedian(abs(eyeAz(preSaccT,apps)),1);
headz=nanmedian(abs(saccAzAll(preSaccT,apps)),1);

eyez=nanmedian(abs(eyeAz(postSaccT,apps)),1);
headz=nanmedian(abs(saccAzAll(postSaccT,apps)),1);

%%% plot errorbar
figure
barweb([headOffset; eyeOffset]',[headOffsetErr; eyeOffsetErr]')
ylabel('azimuth to cricket (deg)');
set(gca,'Xticklabel',{'pre saccade','post saccade'})
legend('head','gaze')
title('Figure 5D: median az pre and post saccade')
%%
clear corr lags corrAll corrAAll uselags 
for vid=1:length(approachEpochs)
    clear hT azC use useN corr corrA 
    nframe = min(length(azimuth{vid}),length(accelerometerChs{vid}(:,6)));
    nframe = min(nframe,length(approachEpochs{vid}));
    useN=approachEpochs{vid}(1:nframe)==0; use=approachEpochs{vid}(1:nframe)==1;
    azC=-(azimuth{vid});
    hT=(accelerometerChs{vid}(:,6))-gyroBias;
    if sum(useN)>3 & sum(~isnan(hT(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corr lags]= nanxcorr(azC(useN),hT(useN),frRate,'coeff');
        hold on;
    else
        corr=NaN;
    end

    if sum(use)>4 & sum(~isnan(hT(use)))>20
        [corrA lagsA]= nanxcorr(azC(use),hT(use),frRate,'coeff');
    else
        corrA=NaN;
    end
    corrAll(vid,:)=corr;
    corrAAll(vid,:)=corrA;
end

figure%('units','normalized','outerposition',[0 0 1 1])
err= nanstd(corrAll)/(sqrt(length(corrAll)));
errA= nanstd(corrAAll)/(sqrt(length(corrAAll)));

shadedErrorBar(1:size(corrAll,2),nanmean(corrAll,1),err,'-b',1); hold on
shadedErrorBar(1:size(corrAAll,2),nanmean(corrAAll,1),errA,'-g',1); hold on
plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 2*frRate+1]); axis square; ylim([-.2 .5]);

clear L
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');
legend(L,{'az to cricket','dtheta head'});

title('Figure 5E: Corr of az and dHead Theta, head is directed to cricket')
%%
clear corr lags corrAll corrAAll uselags mnEye
for vid=1:length(approachEpochs)
    clear hT azC use useN corr corrA
    nframe = min(length(azimuth{vid}),length(accelerometerChs{vid}(:,6)));
    
    nframe = min(nframe,length(approachEpochs{vid}));
    useN=approachEpochs{vid}(1:nframe)==0; use=approachEpochs{vid}(1:nframe)==1;
    
    azC=-(azimuth{vid});
    hT=(accelerometerChs{vid}(:,6))-gyroBias;
    
    mnEye=.5*(thetaR{vid}+thetaL{vid}); mnEye=mnEye-nanmedian(mnEye);
    hT=(accelerometerChs{vid}(:,6))-gyroBias;
    dGazeVid= diff(hT+mnEye);
    
    if sum(useN)>3 & sum(~isnan(dGazeVid(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corr lags]= nanxcorr(azC(useN),dGazeVid(useN),frRate,'coeff');
        hold on;
    else
        corr=NaN;
    end
   
    if sum(use)>4 & sum(~isnan(dGazeVid(use)))>20
        [corrA lagsA]= nanxcorr(azC(use),dGazeVid(use),frRate,'coeff');
    else
        corrA=NaN;
    end
    corrAll(vid,:)=corr;
    corrAAll(vid,:)=corrA;
end

figure%('units','normalized','outerposition',[0 0 1 1])
err= nanstd(corrAll)/(sqrt(length(corrAll)));
errA= nanstd(corrAAll)/(sqrt(length(corrAAll)));

shadedErrorBar(1:size(corrAll,2),nanmean(corrAll,1),err,'-b',1); hold on
shadedErrorBar(1:size(corrAAll,2),nanmean(corrAAll,1),errA,'-g',1); hold on
plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 2*frRate+1]); axis square; ylim([-.2 .5]);

clear L
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');
legend(L,{'az to cricket','dtheta head'});

title('Figure 5F: Corr of az and gaze')
