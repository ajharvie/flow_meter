set(0,'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2);
close all
clear all

% seed = rng, save seed
% load seed
% rng(seed)


%set up signals, filters etc
dt = 0.01;
t = 0:dt:1499*dt;
N = length(t);
fs = 1/dt;
fmin = 1/(t(end)-t(1));
%T = 0.0625;
T = 0.5;
w = 2*pi/T;

%generate signals
phi0 = 0;
% y0 = 0.3*(cos(w*t) + 0.3*(cos(5*(w*t))) + 0.3*(cos(1.4*(w*t))) + 0.2*(cos(0.54*(w*t))) + 0.03*randn(1,N));
% y0 = y0 + 0.55;
% y1 = 0.3*(cos(w*t-pi/3) + 0.3*(cos(5*(w*t-pi/3))) + 0.3*(cos(1.4*(w*t-pi/3))) + 0.2*(cos(0.54*(w*t-pi/3))) + 0.03*randn(1,N));
% y1 = y1 + 0.55;

% y0 = 0.4*(cos(w*t) + 0.00*randn(1,N));
% y0 = y0 + 0.5;
% y1 = 0.4*(cos(w*t-pi/3) + 0.00*randn(1,N));
% y1 = y1 + 0.5;

y0 = smoothts(0.3*square(w*t,75),'g','8') + 0.3*rand(1,N) + 1;
y1 = smoothts(0.3*square(w*t + pi/3,75),'g','8') + 0.3*rand(1,N) + 1;

%digitise signals
nDigi = 12;
y0 = round(y0*(2^nDigi-1))/(2^nDigi-1);
y1 = floor(y1*(2^nDigi-1))/(2^nDigi-1);

%misc parameters
W = hanning(N,'periodic')';
% W = W.^0;
CN = 2;%correction_factor(W,length(W));
k = 0:N-1;
F = linspace(0,fs*(N-1)/N,N);
dF = diff(F); dF = dF(1);
fLo = 0;
fHi = 12;
Fspace = k*dF;

%time domain windowing
Y0a = fft(y0.*W);
Y1a = fft(y1.*W);
phi0a = wrapToPi(angle(Y0a));
phi1a = wrapToPi(angle(Y1a));
R0a = 10*abs(Y0a)/N;
R1a = 10*abs(Y1a)/N;

%find f0
% [~,I] = max(abs(Y0a(3:round(N/2)))); %ignoring DC
% I = I + 2; %correcting for omission of DC
% kEst0a = I - 1 + CN*real((Y0a(I-1)-Y0a(I+1))/(2*Y0a(I)-Y0a(I-1)-Y0a(I+1)));
% [~,I] = max(abs(Y1a(3:round(N/2))));
% I = I + 2;
% kEst1a = I - 1 + CN*real((Y1a(I-1)-Y1a(I+1))/(2*Y1a(I)-Y1a(I-1)-Y1a(I+1)));
% kEsta = (kEst0a+kEst1a)/2;

%find f0 using RMS method
Mboth = sqrt(R0a.*R1a);
[~,I] = max(Mboth(3:round(N/2))); %ignoring DC
I = I + 2;
delta = CN* (Mboth(I+1) - Mboth(I-1))/(2*Mboth(I) + Mboth(I-1) + Mboth(I+1));
kEsta = I - 1 + delta;

%find dphi0a
dphi0a = interp1(1:N,wrapToPi(phi1a-phi0a),kEsta,'linear');

%freq domain windowing
Y0b = sdft_plug(y0);
Y1b = sdft_plug(y1);
phi0b = angle(Y0b);
R0b = 10*abs(Y0b)/N;
phi1b = angle(Y1b);
R1b = 10*abs(Y1b)/N;
dummyA = [NaN Y0b(1:end-1)];
dummyB = [Y0b(2:end) NaN];
Y0c = -0.25*dummyA + 0.5*Y0b - 0.25*dummyB;
R0c = 10*abs(Y0c)/N;
phi0c = angle(Y0c);
dummyA = [NaN Y1b(1:end-1)];
dummyB = [Y1b(2:end) NaN];
Y1c = -0.25*dummyA + 0.5*Y1b - 0.25*dummyB;
R1c = 10*abs(Y1c)/N;
phi1c = angle(Y1c);

%find f0
% [~,I] = max(abs(Y0c(3:round(N/2))));
% I = I + 2;
% kEst0c = I - 1 + CN*real((Y0c(I-1)-Y0c(I+1))/(2*Y0c(I)-Y0c(I-1)-Y0c(I+1)));
% [~,I] = max(abs(Y1c(3:round(N/2))));
% I = I + 2;
% kEst1c = I - 1 + CN*real((Y1c(I-1)-Y1c(I+1))/(2*Y1c(I)-Y1c(I-1)-Y1c(I+1)));
% kEstc = (kEst0c+kEst1c)/2;

%find f0 using RMS method
Mboth = sqrt(R0c.*R1c);
[~,I] = max(Mboth(3:round(N/2))); %ignoring DC
I = I + 2;
delta = CN* (Mboth(I+1) - Mboth(I-1))/(2*Mboth(I) + Mboth(I-1) + Mboth(I+1));
kEstc = I - 1 + delta;

% if abs(kEst1c-kEst0c)/kEst0c > 0.01
%     disp('Help');
% end

%find dphi0c
dphi0c = interp1(1:N,wrapToPi(phi1c-phi0c),kEstc,'linear');


%time domain windowing

figure('units','normalized','outerposition',[0 0 1 1]);
pause(0.5);

axes('position',[0.1 0.85 0.2 0.11]);
plot(t,y0,t,y1);
set(gca,'yLim',[0.5 2]);
xlabel('Time (s)');
yyaxis left
set(gca,'yticklabel',[]);
yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca, 'YColor', 'k');
set(gca,'yLim',yLim);
%set(gca,'yticklabel',[0.5 1]);
text(0.005,1.85,'(a)','fontsize',20);
text(7.5,2.25,'Method 1 - FFT','fontsize',20,'fontweight','bold','horizontalalignment','center');

axes('position',[0.1 0.570 0.2 0.11]);
plot(t,W.*y0,t,W.*y1);
xlabel('Time (s)');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
%yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca, 'YColor', 'k');
%set(gca,'yLim',[0, 0.8]);
%set(gca,'yTick',yTick);
text(0.005,0.9,'(b)','fontsize',20);

axes('position',[0.1 0.30 0.2 0.11]);
plot(Fspace(3:end-1),R0a(3:end-1),'-',Fspace(3:end-1),R1a(3:end-1),'o','markersize',5);
set(gca,'yLim',[0, 0.8]);
line([kEsta*dF kEsta*dF],[0 max(R0a(3:end))],'linestyle','--','color','k','linewidth',1);
set(gca,'xlim',[fLo fHi]);
set(gca,'ylim',[0 0.8]);
set(gca,'xticklabel',[]);
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,0.7,'(c)','fontsize',20);

for i = 1:length(Mboth)
    if Mboth(i)<0.05
        phi0a(i) = 0;
        phi1a(i) = 0;
        phi0b(i) = 0;
        phi1b(i) = 0;
        phi0c(i) = 0;
        phi1c(i) = 0;
    end
end

axes('position',[0.1 0.175 0.2 0.11]);
plot(Fspace(3:end-1),phi0a(3:end-1),'-',Fspace(3:end-1),phi1a(3:end-1),'-','markersize',8);
line([kEsta*dF kEsta*dF],get(gca,'ylim'),'linestyle','--','color','k','linewidth',1);
set(gca,'xlim',[fLo fHi]);
set(gca,'yLim',[-pi, pi]);
%set(gca,'xtick',0:3:18);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
%yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
%set(gca,'yTick',yTick);
text(0.1,2.45,'(d)','fontsize',20);

axes('position',[0.1 0.05 0.2 0.11]);
plot(Fspace(3:end-1),wrapToPi(phi1a(3:end-1)-phi0a(3:end-1)),'-','markersize',8);
line([kEsta*dF kEsta*dF],get(gca,'ylim'),'linestyle','--','color','k','linewidth',1);
xlabel('Frequency (Hz)');
set(gca,'xlim',[fLo fHi]);
set(gca,'yLim',[-pi, pi]);
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\Delta}{\it{\phi}} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.45,'(e)','fontsize',20);

%freq domain filtering

axes('position',[0.375 0.85 0.2 0.11]);
plot(t,y0,t,y1);
set(gca,'yLim',[0.5 2]);
xlabel('Time (s)');
yyaxis left
set(gca,'yticklabel',[]);
yLim = get(gca,'yLim');
ylabel('Signal (a.u.)');
yyaxis right
set(gca, 'YColor', 'k');
set(gca,'yLim',yLim);
%set(gca,'yticklabel',[0.5 1]);
text(0.005,1.85,'(f)','fontsize',20);
text(7.5,2.25,'Method 2 - sDFT','fontsize',20,'fontweight','bold','horizontalalignment','center');

axes('position',[0.375 0.635 0.2 0.11]);
plot(Fspace(3:end-1),R0b(3:end-1),'-',Fspace(3:end-1),R1b(3:end-1),'o','markersize',5);
%set(gca,'xtick',0:3:18);
set(gca,'xticklabel',[]);
set(gca,'xlim',[fLo fHi]);
set(gca,'ylim',[0 1.5]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,1.32,'(g)','fontsize',20);

axes('position',[0.375 0.51 0.2 0.11]);
plot(Fspace(3:end-1),phi0b(3:end-1),'-',Fspace(3:end-1),phi1b(3:end-1),'-','markersize',8);
xlabel('Frequency (Hz)');
set(gca,'xlim',[fLo fHi]);
set(gca,'ylim',[-pi pi]);
%set(gca,'xtick',0:3:18);
%set(gca,'ytick',-1:2:3);
% legend('0','1');
yyaxis left
%set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.45,'(h)','fontsize',20);

axes('position',[0.375 0.30 0.2 0.11]);
plot(Fspace(3:end-1),R0c(3:end-1),'-',Fspace(3:end-1),R1c(3:end-1),'o','markersize',5);
line([kEstc*dF kEstc*dF],[0 max(R0c)],'linestyle','--','Color','k','linewidth',1);
set(gca,'xlim',[fLo fHi]);
set(gca,'ylim',[0 0.8]);
set(gca,'xticklabel',[]);
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\fontname{Times}R} (a.u.)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,0.7,'(i)','fontsize',20);

axes('position',[0.375 0.175 0.2 0.11]);
plot(Fspace(3:end-1),phi0c(3:end-1),'-',Fspace(3:end-1),phi1c(3:end-1),'-','markersize',8);
line([kEstc*dF kEstc*dF],get(gca,'ylim'),'linestyle','--','color','k','linewidth',1);
set(gca,'xlim',[fLo fHi]);
%set(gca,'xtick',0:3:18);
set(gca,'xticklabel',[]);
% legend('0','1');
yyaxis left
set(gca,'yticklabel',[]);
yTick = get(gca,'yTick');
yLim = get(gca,'yLim');
ylabel('{\it\phi} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.45,'(j)','fontsize',20);

axes('position',[0.375 0.05 0.2 0.11]);
plot(Fspace(3:end-1),wrapToPi(phi1c(3:end-1)-phi0c(3:end-1)),'-','markersize',8);
line([kEstc*dF kEstc*dF],get(gca,'ylim'),'linestyle','--','Color','k','linewidth',1);
set(gca,'xlim',[fLo fHi]);
set(gca,'ylim',[-pi pi]);
xlabel('Frequency (Hz)');
yyaxis left
set(gca,'yticklabel',[]);
yLim = get(gca,'yLim');
ylabel('{\Delta}{\it{\phi}} (rad)');
yyaxis right
set(gca,'YColor','k');
set(gca,'yLim',yLim);
set(gca,'yTick',yTick);
text(0.1,2.45,'(k)','fontsize',20);

annotation(gcf,'arrow',[0.125 0.125],[0.837 0.694]);
annotation(gcf,'arrow',[0.125 0.125],[0.837 0.694]-0.275);
annotation(gcf,'arrow',[0.4275 0.4275],[0.835 0.76]);
annotation(gcf,'arrow',[0.4275 0.4275],[0.835 0.76]-0.3375);
annotation(gcf,'textbox',...
    [0.076 0.736 0.048 0.0582],...
    'String',{'apply','window'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'HorizontalAlignment','right');

annotation(gcf,'textbox',...
    [0.095 0.478 0.030 0.032],...
    'String',{'FFT'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.399 0.785 0.0283 0.032],...
    'String',{'sDFT'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right');

annotation(gcf,'textbox',...
        [0.390 0.455 0.035 0.032],...
    'String',{'apply','filter'},...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right');

dummyString = ['$$\hat{f}$$ = ' num2str(kEsta*dF,'%6.5f') ' Hz'];
annotation(gcf,'textbox',...
        [0.185 0.373 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\hat{f}$$ = ' num2str(kEstc*dF,'%6.5f') ' Hz'];
annotation(gcf,'textbox',...
        [0.460 0.373 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\Delta{\hat{\phi}}$$ = ' num2str(dphi0a,'%6.5f')];
annotation(gcf,'textbox',...
        [0.185 0.124 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');

dummyString = ['$$\Delta{\hat{\phi}}$$ = ' num2str(dphi0c,'%6.5f')];
annotation(gcf,'textbox',...
        [0.460 0.124 0.113 0.032],...
    'String',dummyString,...
    'LineStyle','none',...
    'FontSize',20,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'HorizontalAlignment','right',...
    'Interpreter','Latex');


set(gcf,'color','w')

% disp(['k0a = ' num2str(kEst0a,8)]);
% disp(['k1a = ' num2str(kEst1a,8)]);
% disp(['kabar = ' num2str(kEsta,16)]);
% disp(['f0a = ' num2str(dF*kEst0a,8)]);
% disp(['f1a = ' num2str(dF*kEst1a,8)]);
% disp(['fabar = ' num2str(dF*kEsta,8)]);
% disp(['dphia = ' num2str(dphi0a,16)]);
% disp(['% error on dphia = ' num2str(abs(pi/3+dphi0a)/(pi/3))]);
% disp(['Abs error on dphia = ' num2str(pi/3+dphi0a)]);
% disp(['Abs error on dphia in degrees = ' num2str((pi/3+dphi0a)*180/pi)]);

% disp(['k0c = ' num2str(kEst0c,16)]);
% disp(['k1c = ' num2str(kEst1c,8)]);
% disp(['kcbar = ' num2str(kEstc,16)]);
% disp(['f0c = ' num2str(dF*kEst0c,8)]);
% disp(['f1c = ' num2str(dF*kEst1c,8)]);
% disp(['fcbar = ' num2str(dF*kEstc,8)]);
% disp(['dphic = ' num2str(dphi0c,16)]);
% disp(['% error on dphic = ' num2str(abs(pi/3+dphi0c)/(pi/3))]);
% disp(['Abs error on dphic = ' num2str(pi/3+dphi0c)]);
% disp(['Abs error on dphic in degrees = ' num2str((pi/3+dphi0c)*180/pi)]);
