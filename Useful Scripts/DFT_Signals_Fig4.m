set(0, 'defaultAxesFontSize',20);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2);
figure('units','normalized','outerposition',[0.2 0.2 0.4 0.8]);
pause(0.5);

N = length(sig1);
dt = 0.01;
fs = 1/dt;
t = dt*linspace(1,N,N);

W = hanning(N,'periodic');

signal_1 = W.*sig1(1:N);
signal_2 = W.*sig2(1:N);

F = linspace(0,fs*(N-1)/N,N);
dF = diff(F); dF = dF(1);


dft1 = fft(signal_1);
dft2 = fft(signal_2);
mag1 = abs(dft1);
mag2 = abs(dft2);
phase1 = angle(dft1);
phase2 = angle(dft2);


Mboth = sqrt(mag1.*mag2);

[~,I] = max(Mboth(3:round(N/2))); %ignoring DC
%[~,I] = max(Mboth(40:47));
I = I + 2;
delta = 2* (Mboth(I+1) - Mboth(I-1))/(2*Mboth(I) + Mboth(I-1) + Mboth(I+1));
kEst = I - 1 + delta;
fEst = kEst*dF;

phasespec = (phase2-phase1);

for i = 1:N
    if phasespec(i) < 0
        phasespec(i) = phasespec(i) + 2*pi;
    end
    if phase1(i) < 0
        phase1wrapped(i) = phase1(i) + 2*pi;
    end
    if phase2(i) < 0
        phase2wrapped(i) = phase2(i) + 2*pi;
    end
    
end


peakphase = interp1(1:N,wrapToPi(phase2-phase1),kEst,'linear');

if peakphase < 0
    peakphase = peakphase + 2*pi;
end

plotphase = zeros(N,1);
plotphase1 = zeros(N,1);
plotphase2 = zeros(N,1);

for i = 1:N
    if Mboth(i)>150;
        %plotphase(i) = phasespec(i);
        plotphase1(i) = phase1(i);
        plotphase2(i) = phase2(i);
        
%         plotphase1(i) = phase1wrapped(i);
%         plotphase2(i) = phase2wrapped(i);
        plotphase(i) = phasespec(i);
        %plotphase(i) = wrapToPi(phase2(i)-phase1(i));
    end
    
    
end
    
%plotphase = wrapToPi(phase2-phase1);


speed = 2*pi*fEst/peakphase
rate = speed*0.775*60

axpos = -0.46;

axes('position',[0.1 0.85 0.8 0.13]);
plot(t(1:500), sig1(1:500),t(1:500), sig2(1:500))
%plot(t(1:500), sig1(length(sig1)-499:length(sig1)),t(1:500), sig2(length(sig1)-499:length(sig1)))
xlabel('Time (s)');
ylab = ylabel('Signal (a.u.)');
ylab.Position(1) = -0.25;
set(gca,'xLim',[0 5]);
set(gca,'yLim',[1720 1780]);
%set(gca,'xtick',[0 1 2 3 4 5]);
%set(gca,'xticklabel',[0 1 2 13 14 15]);
set(gca,'ytick',[1730 1750 1770]);
set(gca,'yticklabel',[1730 1750 1770]);
rectangle('Position',[2.45,1722,0.1,50],'FaceColor',[1 1 1],'Linestyle','none')
text(2.44,1720,'//','Fontsize',30);
text(2.44,1780,'//','Fontsize',30);
text(0.05, 1775,'(a)','Fontsize',20);
set(gca,'yaxislocation','right');
%ylab.Position(1) = -3;

axes('position',[0.1 0.62 0.8 0.13]);
plot(F(3:length(F)), mag1(3:length(F)), F(3:length(F)), mag2(3:length(F)))
set(gca,'yLim',[0 2000]);
set(gca,'xLim',[0 10]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[0 2000]);
set(gca,'yticklabel',[0 2000]);
ylab = ylabel('$$R$$ (a.u.)', 'Interpreter', 'LaTeX');
ylab.Position(1) = axpos;
text(0.093, 1810,'(b)','Fontsize',20);
set(gca,'yaxislocation','right');

axes('position',[0.1 0.47 0.8 0.13]);
plot(F(3:length(F)), plotphase1(3:length(F)), F(3:length(F)), plotphase2(3:length(F)))
%plot(F, phasespec)
set(gca,'yLim',[-2*pi 2*pi]);
set(gca,'yTick',[-2*pi, 0, 2*pi]);
set(gca,'yTickLabel',{'-2\pi', '0', '2\pi'});
set(gca,'xLim',[0 10]);
set(gca,'xticklabel',[]);
ylab = ylabel('$$\tilde{\phi}$$', 'Interpreter', 'LaTeX');
ylab.Position(1) = axpos;
text(0.093, 5.089,'(c)','Fontsize',20);
set(gca,'yaxislocation','right');

axes('position',[0.1 0.32 0.8 0.13]);
% plot(F, phase1, F, phase2)
plot(F(3:length(F)), plotphase(3:length(F)))
set(gca,'yLim',[0 2*pi]);
set(gca,'yTick',[0, 2*pi]);
set(gca,'yTickLabel',{'0', '2\pi'});
set(gca,'xLim',[0 10]);
xlabel('Frequency (Hz)');
ylab = ylabel('$$\Delta\tilde{\phi}$$', 'Interpreter', 'LaTeX');
ylab.Position(1) = axpos;
xline(fEst,'--','LineWidth',2);
text(7.2,4.9,'$$\Delta\tilde{\phi}$$ = 3.1217 rad', 'Interpreter', 'LaTeX','fontsize',20);
text(0.093, 5.68,'(d)','Fontsize',20);
set(gca,'yaxislocation','right');


axes('position',[0.1 0.08 0.8 0.13]);
plot(F(3:length(F)), Mboth(3:length(F)))
set(gca,'yLim',[0 1300]);
set(gca,'xLim',[0 10]);
%set(gca,'xticklabel',[]);
xlabel('Frequency (Hz)');
ylab = ylabel('$$\tilde{R}_{AB}$$ (a.u.)', 'Interpreter', 'LaTeX','fontsize',20);
ylab.Position(1) = axpos;
xline(fEst,'--','LineWidth',2);
text(0.093, 1190,'(e)','Fontsize',20);
set(gca,'yaxislocation','right');
text(7.5,1050,'$$\hat{f}$$ = 2.1311 Hz', 'Interpreter', 'LaTeX','fontsize',20);



% subplot(5,1,1)
% hold on
% plot(t, sig1(1:N), t, sig2(1:N))
% subplot(5,1,2)
% plot(F,mag1)
% subplot(5,1,3)
% plot(F,mag2)
% subplot(5,1,4)
% plot(F,Mboth)


