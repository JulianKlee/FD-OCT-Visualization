clear all; close all;format long;

IRabsorb=[954.3, 3.061;
956.9, 3.173;
960.8, 3.286;
964.7, 3.411;
969.9, 3.486;
978.9, 3.475;
984.0, 3.401;
989.0, 3.252;
992.8, 3.153
997.8, 2.980;
1003, 2.819;
1008, 2.645;
1012, 2.459;
1017, 2.273;
1020, 2.087;
1024, 1.914;
1028, 1.740;
1034, 1.567;
1039, 1.381;
1047, 1.233;
1056, 1.122;
1066, 1.061;
1076, 1.087;
1084, 1.163;
1090, 1.263;
1096, 1.376;
1099, 1.476;
1102, 1.576;
1105, 1.700;
1108, 1.875;
1110, 1.999;
1113, 2.298;
1116, 2.684;
1120, 3.156;
1127, 3.965;
1132, 4.662;
1146, 6.081;
1157, 7.201;
1172, 8.856;
1180, 9.516];
%lambeer=4pi*k(lambda)/lambda ; values are ke-6 und lambda is 1e-6 kürzt
%sich also

%% x vectors
xLambda=linspace(0, 2000,2000);
xNormalized=linspace(0,1,2000);

%% spectral shaping of the eye

eyeSpecShaper = 1./(interp1(IRabsorb(:,1),IRabsorb(:,2),xLambda,'spline'));
eyeSpecShaper(1:950)=0;eyeSpecShaper(1200:end)=0;

eyeSpecShapeNaN=eyeSpecShaper;
eyeSpecShapeNaN(eyeSpecShapeNaN==0)=nan;

figure;
plot(xLambda,eyeSpecShapeNaN,'color',[0.4660,0.6740,0.1880]','LineWidth',2);
title('interpol eye spec shape');
xlim([950 1200]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'interpol eye spec shape','svg');


%% rectangular spectra

bigSpectrumRect=zeros(1,2000);
bigSpectrumRect(1000:1150)=1;
rectangleBig=bigSpectrumRect;

smallSpectrumRect=zeros(1,2000);
smallSpectrumRect(1050:1100)=1;
rectangleSmall=smallSpectrumRect;

figure;
h=plot(xLambda,rectangleBig,xLambda,rectangleSmall,'--','LineWidth',2);
title('rectangle spectra');
set(h,{'color'},{[0.8500, 0.3250, 0.0980];	[0, 0.4470, 0.7410]});
xlim([950 1200]);ylim([-0.1 1.1]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'rectangle spectra','svg');

%% rectangular spectra formed by the eye

bigSpectrumEye=bigSpectrumRect.*eyeSpecShaper;
bigSpectrumEyeNaN=bigSpectrumEye;
bigSpectrumEyeNaN(bigSpectrumEyeNaN==0)=nan;

smallSpectrumEye=smallSpectrumRect.*eyeSpecShaper;
smallspectrumEyeNaN=smallSpectrumEye;
smallspectrumEyeNaN(smallspectrumEyeNaN==0)=nan;

figure;
h=plot(xLambda,bigSpectrumEyeNaN,xLambda,smallspectrumEyeNaN,'--','LineWidth',2);
set(h,{'color'},{[0.8500, 0.3250, 0.0980];	[0, 0.4470, 0.7410]});
title('shaped spetra');
xlim([950 1200]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'shaped spectra','svg');

%% ifft of the rect spectra

ifftBigSpectrumRect=abs(ifft(rectangleBig.^2));%square acording to formula diss W.W.
ifftBigSpectrumRect=1/max(ifftBigSpectrumRect)*ifftBigSpectrumRect;

% figure;
% plot(xNormalized,abs(ifftBigSpectrumRect));
% title('ifft abs of the big Recta spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

ifftSmallSpectrumRect=abs(ifft(smallSpectrumRect.^2));
ifftSmallSpectrumRect=1/max(ifftSmallSpectrumRect)*ifftSmallSpectrumRect;

% figure;
% plot(xNormalized,ifftSmallSpectrumRect);
% title('ifft abs of small Rect spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

figure;
plot(xNormalized,ifftSmallSpectrumRect,xNormalized,ifftBigSpectrumRect,'LineWidth',2);
title('ifft abs of the two Rect spectra zoom');
xlim([0 0.25]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'ifft abs of the two Rect spectra zoom','svg');

figure;
plot(xNormalized,ifftSmallSpectrumRect,xNormalized,ifftBigSpectrumRect,'LineWidth',2);
title('ifft abs of the two Rect spectra');
xlim([0 1]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'ifft abs of the two Rect spectra','svg');

%% fft of the rect spectra

fftBigSpectrumRect=abs(fft(bigSpectrumRect.^2));%square acording to formula diss W.W.
fftBigSpectrumRect=1/max(fftBigSpectrumRect)*fftBigSpectrumRect;

% figure;
% plot(xNormalized,fftBigSpectrumRect);
% title('fft abs of the big Rect spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

fftSmallSpectrumRect=abs(fft(smallSpectrumRect.^2));
fftSmallSpectrumRect=1/max(fftSmallSpectrumRect)*fftSmallSpectrumRect;

% figure;
% plot(xNormalized,fftSmallSpectrumRect);
% title('fft abs of small Rect spectrum');
% xlim([0 1]);pbaspect([1 1 1]);
% 
% figure;
% plot(xNormalized,fftSmallSpectrumRect,xNormalized,fftBigSpectrumRect);
% title('fft abs of the two Rect spectra');
% xlim([0 0.25]);pbaspect([1 1 1]);

%% ifft of the eye shaped spectra

ifftBigSpectrumEye=abs(ifft(bigSpectrumEye.^2));%square acording to formula diss W.W.
ifftBigSpectrumEye=1/max(ifftBigSpectrumEye)*ifftBigSpectrumEye;

% figure;
% plot(xNormalized,abs(ifftBigSpectrumEye));
% title('ifft abs of the big spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

ifftSmallSpectrumEye=abs(ifft(smallSpectrumEye.^2));
ifftSmallSpectrumEye=1/max(ifftSmallSpectrumEye)*ifftSmallSpectrumEye;
% 
% figure;
% plot(xNormalized,ifftSmallSpectrumEye);
% title('ifft abs of small spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

figure;
plot(xNormalized,ifftSmallSpectrumEye,xNormalized,ifftBigSpectrumEye,'LineWidth',2);
title('ifft abs of the two spectra');
xlim([0 0.25]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'ifft abs of the two eye spectra zoom','svg');

figure;
plot(xNormalized,ifftSmallSpectrumEye,xNormalized,ifftBigSpectrumEye,'LineWidth',2);
title('ifft abs of the two spectra');
xlim([0 1]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'ifft abs of the two eye spectra','svg');
%% fft of the eye shaped spectra

fftBigSpectrumEye=abs(fft(bigSpectrumEye.^2));%square acording to formula diss W.W.
fftBigSpectrumEye=1/max(fftBigSpectrumEye)*fftBigSpectrumEye;

% figure;
% plot(xNormalized,fftBigSpectrumEye);
% title('fft abs of the big spectrum');
% xlim([0 1]);pbaspect([1 1 1]);

fftSmallSpectrumEye=abs(fft(smallSpectrumEye.^2));
fftSmallSpectrumEye=1/max(fftSmallSpectrumEye)*fftSmallSpectrumEye;

% figure;
% plot(xNormalized,fftSmallSpectrumEye);
% title('fft abs of small spectrum');
% xlim([0 1]);pbaspect([1 1 1]);
% 
% figure;
% plot(xNormalized,fftSmallSpectrumEye,xNormalized,fftBigSpectrumEye);
% title('fft abs of the two spectra');
% xlim([0 0.25]);pbaspect([1 1 1]);

%% modulation sin
frq=100;
for ki=1:3
    moduSin(ki,:)=sin(frq*2*pi*xNormalized);

    figure;
    plot(xNormalized,moduSin(ki,:),'color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
    title(strcat('modulation Sin',num2str(ki)));
    xlim([0 0.15]);ylim([-1.1 1.1]);pbaspect([1 1 1]);
    frq=frq+100;
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('modulation Sin',num2str(ki)),'svg');
end
%% ifft moduSin

for ki=1:3
    dirac(ki,:)=abs(ifft(moduSin(ki,:)));

    figure;
    plot(xNormalized,dirac(ki,:),'color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
    title(strcat('dirac',num2str(ki)));
    pbaspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('dirac',num2str(ki)),'svg');
end

%% convolve transformed rect spectrum and moduSin
for ki=1:3
    ConvBigSpecShapeRect=conv(ifftBigSpectrumRect,dirac(ki,:));
    ConvBigSpecShapeRect=ConvBigSpecShapeRect(1:2000);%remove extra from convolution
    ConvBigSpecShapeRect=1/max(ConvBigSpecShapeRect)*ConvBigSpecShapeRect;

    ConvSmallSpecShapeRect=conv(ifftSmallSpectrumRect,dirac(ki,:));
    ConvSmallSpecShapeRect=ConvSmallSpecShapeRect(1:2000);%remove extra from convolution
    ConvSmallSpecShapeRect=1/max(ConvSmallSpecShapeRect)*ConvSmallSpecShapeRect;

    figure;
    plot(xNormalized,ConvSmallSpecShapeRect,xNormalized, ConvBigSpecShapeRect,'LineWidth',2);
    title(strcat('convolve transformed rect spectrum and moduSin ',num2str(ki)));
    xlim([0 0.25]);pbaspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('convolve transformed rect spectrum and moduSin zoom',num2str(ki)),'svg');
    
    figure;
    plot(xNormalized,ConvSmallSpecShapeRect,xNormalized, ConvBigSpecShapeRect,'LineWidth',2);
    title(strcat('convolve transformed rect spectrum and moduSin ',num2str(ki)));
    pbaspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('convolve transformed rect spectrum and moduSin',num2str(ki)),'svg');
end
%% convolve transformed eye spectrum and moduSin
for ki=1:3
    ConvBigSpecShapeEye=conv(ifftBigSpectrumEye,dirac(ki,:));
    ConvBigSpecShapeEye=ConvBigSpecShapeEye(1:2000);%remove extra from convolution
    ConvBigSpecShapeEye=1/max(ConvBigSpecShapeEye)*ConvBigSpecShapeEye;
    
    ConvSmallSpecShapeEye=conv(ifftSmallSpectrumEye,dirac(ki,:));
    ConvSmallSpecShapeEye=ConvSmallSpecShapeEye(1:2000);%remove extra from convolution
    ConvSmallSpecShapeEye=1/(max(ConvSmallSpecShapeEye))*ConvSmallSpecShapeEye;

    figure;
    plot(xNormalized,ConvSmallSpecShapeEye,xNormalized, ConvBigSpecShapeEye,'LineWidth',2);
    title(strcat('convolve transformed eye spectrum and moduSin zoom',num2str(ki)));
    xlim([0 0.25]);pbaspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('convolve transformed eye spectrum and moduSin zoom',num2str(ki)),'svg');
    
    figure;
    plot(xNormalized,ConvSmallSpecShapeEye,xNormalized, ConvBigSpecShapeEye,'LineWidth',2);
    title(strcat('convolve transformed eye spectrum and moduSin ',num2str(ki)));
    pbaspect([1 1 1]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
    saveas(gcf,strcat('convolve transformed eye spectrum and moduSin',num2str(ki)),'svg');
    

end
%% convolve with fourier space multiplication: rect

FFTBigSpecShapeRect=abs(fft(moduSin(2,:).*rectangleBig));
FFTBigSpecShapeRect=1/max(FFTBigSpecShapeRect)*FFTBigSpecShapeRect;

FFTSmallSpecShapeRect=abs(fft(moduSin(2,:).*rectangleSmall));
FFTSmallSpecShapeRect=1/max(FFTSmallSpecShapeRect)*FFTSmallSpecShapeRect;

figure;
plot(xNormalized,moduSin(2,:).*rectangleBig,'color',[0.8500, 0.3250, 0.0980]','LineWidth',2);
title('multiply rect big * moduSin');
xlim([0.4 0.7]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply rect big moduSin','svg');

figure;
plot(xNormalized,moduSin(2,:).*rectangleSmall,'LineWidth',2);
title('multiply rect small * moduSin');
xlim([0.4 0.7]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply rect small moduSin','svg');

figure;
plot(xNormalized,FFTSmallSpecShapeRect,xNormalized,FFTBigSpecShapeRect,'LineWidth',2);
title('multiply in F: rect big * moduSin');
pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply in F rect moduSin','svg');

figure;
plot(xNormalized,FFTSmallSpecShapeRect,xNormalized,FFTBigSpecShapeRect,'LineWidth',2);
title('multiply in F: rect big * moduSin');
xlim([0 0.25]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply in F rect moduSin zoom','svg');

%% convolve with fourier space multiplication: eye shaped

FFTeyeBigSpecShapeRect=abs(fft(moduSin(2,:).*bigSpectrumEye));
FFTeyeBigSpecShapeRect=1/max(FFTeyeBigSpecShapeRect)*FFTeyeBigSpecShapeRect;

FFTeyeSmallSpecShapeRect=abs(fft(moduSin(2,:).*smallSpectrumEye));
FFTeyeSmallSpecShapeRect=1/max(FFTeyeSmallSpecShapeRect)*FFTeyeSmallSpecShapeRect;

figure;
plot(xNormalized,moduSin(2,:).*bigSpectrumEye,'color',[0.8500, 0.3250, 0.0980]','LineWidth',2);
title('multiply eye big * moduSin');
xlim([0.4 0.7]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply eye big moduSin','svg');

figure;
plot(xNormalized,moduSin(2,:).*smallSpectrumEye,'LineWidth',2);
title('multiply eye small * moduSin');
xlim([0.4 0.7]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply eye small moduSin','svg');

figure;
plot(xNormalized,FFTeyeSmallSpecShapeRect,xNormalized,FFTeyeBigSpecShapeRect,'LineWidth',2);
title('multiply in F: eye rect Big * moduSin');
pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply in F eye rect Big moduSin','svg');

figure;
plot(xNormalized,FFTeyeSmallSpecShapeRect,xNormalized,FFTeyeBigSpecShapeRect,'LineWidth',2);
title('multiply in F: eye rect Big * moduSin');
xlim([0 0.25]);pbaspect([1 1 1]);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'Visible','off');
saveas(gcf,'multiply in F eye rect Big moduSin zoom','svg');
%%



