%simulator
t_delay = 0.15; %delay time in seconds
f = 3;      %frequency of oscillations

duty = 80;  % duty cycle of wave 
timestep = 0.01;

ph = 2*pi*f*t_delay;

t = timestep*linspace(0, 99999,100000);
s1 = 1*smoothdata(0.8*square(f*(2*pi*t), duty) + 1.5,'gaussian') + 0.01*randn(1,100000);
s2 = 1*smoothdata(0.8*square(f*(2*pi*t + ph/f), duty) + 1.5,'gaussian') + 0.01*randn(1,100000);


y1 = s1;
y2 = s2;

plot(t(1:500), s1(1:500),t(1:500), s2(1:500))


