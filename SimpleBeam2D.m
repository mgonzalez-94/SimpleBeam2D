% SimpleBeam2D
% 
% Response a simple supported beam in 2D.
% 
% INPUTS:
% 
% L: length.
% E: Young's modulus.
% g: mass per unit length.
% Iz: inertia around z axis (m4).
% ndof: number of degrees of freedom.
% fs: sample frequency (Hz)
% Amp: amplitude of the inv(M)*F excitation
% tf: time length (s)
% 
% OUTPUTS:
% 
% data: text file with [t,Out.y] (s) (m/s2)
% 
% Units: Internaional system.
% 
% %%%%%%%%%%%%%%%%%%%
% %%% Mateo G. H. %%%
% %%% 2021/05/09  %%%
% %%%%%%%%%%%%%%%%%%%
clc, clearvars, close all;
%% Structure
%%% -----------------------------------------------------------------------
% Mechanical properties
St.L  = 60; % (m)
St.E  = 25*1e9; % (Pa)
St.g  = 1529; % (kg/m)
St.Iz = 6.75; %(m4)
St.ndof = 6; % number of dof
St.n  = (1:St.ndof)';
St.z  = 0.05*ones(St.ndof,1); % Damping ratio
St.Z  = diag(St.z); % Damping ratio
%%% -----------------------------------------------------------------------
St.x  = linspace(0,St.L,St.ndof+2); % x divisions
St.x([1,end]) = [];
St.Ph  = sin(St.n.*pi.*St.x./St.L);
St.iPh = inv(St.Ph);
St.wn  = ((St.n.*pi./St.L).^2).*sqrt(St.E*St.Iz/St.g);
St.Wn  = diag(St.wn);
St.Wn2 = St.Wn.^2;
St.fn  = St.wn./(2*pi);
%%% -----------------------------------------------------------------------
St.Ass = [zeros(St.ndof),eye(St.ndof);-St.Wn2, -2*St.Z*St.Wn];
St.Bss = [zeros(St.ndof);St.iPh];
St.Css = [-St.Ph*St.Wn2, -St.Ph*2*St.Z*St.Wn];
St.Dss = 0;
St.Sys = ss(St.Ass,St.Bss,St.Css,St.Dss);
%% Input
In.fs   = 500;              % (Hz)
In.amps = 0.01*9.81;        % (m/s2)
In.dt   = 1/In.fs;          % (s)
In.tf   = 60*1;            % (s)
In.t    = (0:In.dt:In.tf)'; % (s)
In.nsam = length(In.t);     % Número de muestras
%%% Noise
rng('default');
In.F = In.amps*rand(In.nsam,St.ndof)-In.amps/2; % (m/s2)
In.F([1,end],:) = 0;
In.ncha = size(In.F,2);
%%% Frecuencia
In.L = length(In.F);                       % Longitud de la señal
In.Y = 2*(fft(In.F)./In.L);                % (Y(f)) transformada rápida de Fourier.
In.f = ((0:length(In.Y)/2-1)*In.fs/In.L)'; % (Hz)
In.Y = In.Y(1:length(In.f),:);             % (|Y(f)|)
%% Out
%%% -----------------------------------------------------------------------
% Respuesta sin ruido, aceleración total en coord. geom.
Out.y = lsim(St.Sys,In.F,In.t);  % (m/s2)
%%%
data{1,1} = [In.t,Out.y];
%%% Add noise to signal ---------------------------------------------------
rng('default');
Out.snr = 5:40; % (dB)
Lsnr = length(Out.snr);
Out.yn = cell(Lsnr,1);
for ii = 1:Lsnr
    % Respuesta con ruido adicionado
    Out.yn{ii} = awgn(Out.y,Out.snr(ii)); % (m/s2)
    data{1+ii,1} = [In.t,Out.yn{ii}];
end
%% Export
%%% -----------------------------------------------------------------------
FileName = num2cell((1:Lsnr+1)');
FileName = cellfun(@(Cell) char(['data_',num2str(Cell,'%.0f'),'.txt']),FileName,'uni',false);
%%% -----------------------------------------------------------------------
NumberFormat = char(ones(size(data{1},2),1)*' %10.4f')';
FileFormat   = [NumberFormat(:)','\r\n'];
for ii=1:length(data)
    FileID = fopen(FileName{ii},'w');
    fprintf(FileID,FileFormat,data{ii}');
    fclose(FileID);
end