clear all;close all; format long;

filename='ssoct1.gif';
filename2='ssoct';

amplitude=.1;
samplingFrq=1/1e5;
frqStart=10;
frqEnd=40;
phaseOffset=90;%degrees;
t0=0;
tSweepDuration=2;

t1=t0+tSweepDuration;

tChirp = t0:samplingFrq:t1-samplingFrq;
tSweepDurationLength=length(tChirp);
yChirp = amplitude*chirp(tChirp,frqStart,2,frqEnd,'linear',phaseOffset);

% figure;
% plot(tChirp,yChirp);
% ylim([-2 2]); xlim([-0.5 2.5]);

tTotalDuration=4;
tOffsetStart=1;
tOffsetEnd=1;

tTotal=t0:samplingFrq:tTotalDuration-samplingFrq;
LtTotal=length(tTotal);

offsetSecChirp=2;% from 0 to 2

yPrimarySweep=nan(1,LtTotal);
yPrimarySweep(1/4*LtTotal+1:3/4*LtTotal)=yChirp;

figure;
plot(tTotal,yPrimarySweep);
ylim([-2 2]); xlim([-0.5 4.5]);

figure;
offsetVec=0.6:0.001:1;
tTot4tel=linspace(0,1,.25*LtTotal);
LoffsetVec=length(offsetVec);
%amplitudeVec=(abs(sin(linspace(0,2*pi,length(offsetVec))))+0.2)*(1/1.2);
amplitudeVec=zeros(1,length(offsetVec));
amplitudeVec(round(length(offsetVec)*.1))=1;
amplitudeVec(round(length(offsetVec)*.25))=.75;
amplitudeVec(round(length(offsetVec)*.32))=1;
amplitudeVec(round(length(offsetVec)*.5))=.3;
amplitudeVec(round(length(offsetVec)*.75))=1;
amplitudeVec(round(length(offsetVec)*.95))=.15;
figure;plot(amplitudeVec);
amplitudeVec=conv(conv(amplitudeVec,gausswin(100),'same'),gausswin(50),'same');
figure;plot(amplitudeVec);
ampliVecT4tel=interp1(linspace(0,1,length(amplitudeVec)),amplitudeVec,tTot4tel,'pchip');
figure;plot(ampliVecT4tel);

% ld=length(amplitudeVec)-LoffsetVec;
% amplitudeVec=amplitudeVec(round(ld/2):round(ld/2)+LoffsetVec);
fftEnv=zeros(length(amplitudeVec),.25*LtTotal);
figure('units','normalized','outerposition',[0 0 1 1]);
for k=1:length(amplitudeVec)
    tmp1=zeros(1,length(amplitudeVec));
    tmp1(k)=amplitudeVec(k);
    tmp1=interp1(linspace(0,1,length(amplitudeVec)),tmp1,tTot4tel);
    tmp1=conv(tmp1,gausswin(5000),'same');
%     plot(tmp1);hold on;pause(.5);
    fftEnv(k,:)=tmp1;
end

fftEnv=1/max(fftEnv(:))*fftEnv*.5;
ampliVecT4tel=1/max(ampliVecT4tel(:))*ampliVecT4tel*.5;
amplitudeVec=1/max(amplitudeVec)*amplitudeVec;

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(offsetVec)
    
    offsetSecChirp=offsetVec(i);
    amplitudeSecChirp=amplitudeVec(i);
    
    ySecondarySweep=nan(1,LtTotal);
    ySecondarySweep(round(offsetSecChirp*1/4*LtTotal)+1:round(offsetSecChirp*1/4*LtTotal)+tSweepDurationLength)=amplitudeSecChirp*yChirp;
    
    frqOverlap=0.5*(ySecondarySweep+yPrimarySweep);
    frqOverlapsqrd=(frqOverlap).^2;
    
    NotIsnanIndices=find(~isnan(frqOverlap));
    firstElement=min(NotIsnanIndices);
    lastElement=max(NotIsnanIndices);
    lengthTEnvelope=length(NotIsnanIndices);
    tEnvelope=linspace(tTotal(firstElement),tTotal(lastElement),lengthTEnvelope);
    
    frqOverlapsrd0=frqOverlapsqrd;
    frqOverlapsrd0(isnan(frqOverlapsrd0))=[];
    [upEnvsqrd,loEnvsqrd]=envelope(frqOverlapsrd0,3000,'peak');
    
    a=abs(fft(upEnvsqrd));
    [w,ind]=max(fftEnv(i,:));

    
    frqOverlap0=frqOverlap;
    frqOverlap0(isnan(frqOverlap0))=[];
    [upEnv,loEnv]=envelope(frqOverlap0,300,'peak');
    upEnv(upEnv>.1)=nan;upEnv(upEnv<0)=nan;
    loEnv(loEnv<-.1)=nan;loEnv(loEnv>0)=nan;
    
    h=plot(tTotal,yPrimarySweep+0.3,tTotal,ySecondarySweep,...
        tTotal+2.5,frqOverlap+0.175,tEnvelope+2.5,upEnv+.175,tEnvelope+2.5,loEnv+.175,...
        tTotal+2.5,50*frqOverlapsqrd-.75,tEnvelope+2.5,50*upEnvsqrd-.75,...
        tTot4tel+1.5,fftEnv(i,:)-0.75,tTot4tel+1.5,ampliVecT4tel-0.75,...
        tTot4tel(ind)+1.5,w-0.75,'o','Linewidth',2);
    set(h,{'color'},{
        [0, 0.4470, 0.7410]; [0.9290, 0.6940, 0.1250];
        [0.4660, 0.6740, 0.1880];[0.8500, 0.3250, 0.0980];
        [0.8500, 0.3250, 0.0980];[0.4660, 0.6740, 0.1880]*.7;
        [0.8500, 0.3250, 0.0980];[0.4940, 0.1840, 0.5560];
        [0.6350, 0.0780, 0.1840];[0.6350, 0.0780, 0.1840]});
    h=gca;
    set(h,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    xlim([.5 5.5]);ylim([-0.8 .5]);
    pbaspect([2 1 1]);
    drawnow; 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      saveas(gcf,strcat(filename2,num2str(i)),'svg');
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',0.03, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.03,'WriteMode','append'); 

      end 
    

    
end



