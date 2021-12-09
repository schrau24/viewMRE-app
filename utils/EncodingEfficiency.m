function enceff=EncodingEfficiency( f, fMEG, G, dt, M1, time )
%% EncodingEfficiency(f,fMEG,G,dt) of a Bipolar MEG
% Calculates the encoding efficiency of a bipolar meg with finite slew rate.
% 
% Inputs:
%   f       wave frequency  [Hz]
%   fMEG    MEG frequency   [Hz]
%           = 1 / duration of the bipolar motion encoding gradient
%   G       gradient strength [mT/m]
%   dt      duration to slew to maximum gradient strength [ms] ( = G / slew rate )
%   M1      if true, calculates encoding efficiency for tripolar M1 nulled
%           MEG. otherwise assumes bipolar MEG
%   (time)  optional: 
%           the time at which the MEG is played is defined 
%           relative to an arbitrary time point 
%           (usually I chose the RF pulse time, and the mid of the bipolar gradient)
%           This is only necessary if wave fields are to be compared between sequences 
%           with different encoding schemes (e.g. Unipolar MRE) or different MEG placement
%           with respect to the trigger of the wave generator.
%
% Output:   Complex encoding efficiency in rad / micro meter
% 
% References:
% Guenthner, C. et al. Unipolar MR elastography: Theory, numerical analysis and implementation. NMR Biomed. 33, e4138 (2020).
% Guenthner, C. & Kozerke, S. Encoding and readout strategies in magnetic resonance elastography. NMR Biomed. 31, e3919 (2018).
%
% Christian Guenthner, 2021 - guenthner@biomed.ee.ethz.ch
% Eric Schrauben, 2021 - e.m.schrauben@amsterdamumc.nl

% time: set to zero if not provided
if nargin<6 || isempty(time); time = 0; end
% gyromagnetic ratio
gamma = 2*pi*42.576;

% convert frequencies to kHz to be consistent with ms
f    = f / 1000; %[kHz]
fMEG = fMEG / 1000; %[kHz]
omega = f * 2 * pi;

% duration of the MEG
T = 1./ fMEG; %[ms]

% result using int G exp(i omega t) definition
if M1
    % for tripolar (M1 = 0 MEG), eric 20210208, taken from Dittman et al.
    % 2018, In Vivo Multifrequency wMRE of the Human Brain and Liver, Eqn 5
    q = T*f; % encoding fraction (MEG period * transducer frequency)
%     sinfacs = sin(pi*q) - 2*sin(3/2*pi*q) + sin(2*pi*q);
%     enceff = gamma * G * (1/f) * (1/pi) * sinfacs;
    
    % using the calculation from Guenther et al.
    beta = f*dt;
    sinfacs = cos(2*pi*beta) + sin(3/2*pi*beta) + sin(5/2*pi*beta);
    enceff = (gamma*abs(sinfacs))/(2*pi^2*f*beta) * G;
else
    % for bipolar (M0 = 0 MEG)
    sinfacs = sin(T*omega/4) .* sin(dt*omega/2) .* sin(dt*omega/2-T*omega/4);
    enceff = 8*1i * G./(dt .* omega.^2) .* gamma .* exp(1i*(T.*omega)/2)*sinfacs;
end
% correct for time shift of meg
enceff = enceff .* exp(1i*time*omega); %
enceff = enceff*1e-6; %to [rad/um]

% print result
fprintf('Enc. Efficiency = %f rad/um\n',enceff);