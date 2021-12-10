% makeSARdata
%  make SAR project data
%  Written by Sajid Ahmed
clear all,  clc
c = 3e8;                               % velocity of light in m/sec
%% ---------------------------------------------------
% User input section  
file = 'MyFile.dat';
Dcr = 10;           % cross-range resolution, m
DR = 10;            % range resolution, m
Ro = 30000;         % range to swath center Slant Range
v = 150;            % platform velocity, m/s (540km/h)
fc = 10e9;          % RF frequency in Hz
lambda = c/fc;      % wavelength, m
T = 5e-6;           % Chirp pulse duration, 
Daz = 0.2;          % antenna azimuth length, m
Qaz = lambda/Daz;   % azimuth beamwidth, radians
Ls = 3000;          % swath depth, m
BWdopp = 2*v/lambda*Qaz;     % slow time Doppler bandwidth, Hz or Maximum Doppler Shift
PRF = 5*BWdopp; PRI = 1/PRF;
B = c/(2*DR);
Fs = 3*B; Ts = 1/Fs;
tf = 0:Ts:T;
st = exp(j*pi*(B/T)*(tf-T/2).^2);  % Transmitted Signal

Ta = Ro*lambda/(2*v*Dcr);
Dsar = v*Ta;
Qsar = lambda/(2*Dsar);

Nc = round(Ta/PRI);
ts = linspace(-Ta/2,Ta/2,Nc);    % Start of slow each pulse time     
u  = linspace(-(Dsar/2),(Dsar/2),Nc);  % Distance of SAR from the reference point on the direction of motion path.  

% Ranges are relative to the CRP range (Rp) above.
% coords = [0,0];         % a single point target at the CRP
TarLocs = ...              % a grid of 9 point targets
        [-1000,-1000;
         -1000,0;
         -1000,+1000;
             0,-1000;
             0,0;
             0,+1000;
         +1000,-1000;
         +1000,0;
         +1000,+1000];
%% Find the range of SAR at each location, u, along the flight from each target.
for tt = 1:length(TarLocs)
    R(tt,:) = sqrt((Ro-TarLocs(tt,2))^2 + (TarLocs(tt,1) - u).^2);
end
% At this stage the range of Tth target at each location u is in mth row. 

%% Now we have to find the received signal at fast time samples 
% For swath length Ls, the min and max ranges will be (Ro-Ls/2) and (Ro+Ls/2) 
%% Range of fast times
tf = 2*(Ro-Ls/2)/c:Ts:(2*(Ro+Ls/2)/c + T); 
ind_tf = int32(tf*Fs);
tmin = 2*(Ro-Ls/2)/c;
Y = zeros(length(tf),Nc);

for p = 1:Nc
    for tt = 1:length(TarLocs)
        tau = 2*R(tt,p)/c;
        t = tau:Ts:(tau+T);
        tind = int32(Fs*(t-tmin));
        Y(tind,p) = Y(tind,p) + exp(-j*2*pi*fc*tau)*exp(j*pi*B/T*(t'-T/2-tau).^2);
        %% pth col of Y contains received data from different targets after pth pulse 
    end
end
figure(1)
imagesc(1:Nc,tf,real(Y));
xlabel('Pulse Number (p)')
ylabel('Time')
title('Received signal at different different time after pulse p')
%%        
%% Processing at the Receiver
h = conj(fliplr(st));
for p = 1:Nc
    Ym(:,p) = conv(h,Y(:,p));
end
imagesc(abs(Ym))

keyboard





     