%% Orenthal Server
addpath('/nfs/37zfs1-packages/Matlab/packages') % QUIUC and BIGQUIC path

%% on my Mac
addpath(genpath('/Users/sinajafarzadeh/Desktop/MLCB_extension/QUIC')) % QUIC ...
addpath(genpath('/Users/sinajafarzadeh/Desktop/MLCB_extension/mlcb')) % Project Path

%% Windows at home
addpath(genpath('C:\Users\jafar\Desktop\L1precision')) % QUIC library in our machine
addpath(genpath('C:\Users\jafar\Desktop\KoborDNAm'))
addpath(genpath('../../KoborDNAm')) % dataset folder
addpath(genpath('../../../KoborDNAm')) % dataset folder
addpath(genpath('../../../UCSC')) % dataset folder



% date_time_now_str=datestr(datetime('now'));
% date_time_now_str = strrep(date_time_now_str,':','-');
% logger_file_name=strcat(mfilename,'_',date_time_now_str,'.txt');
% logger = log4m.getLogger(logger_file_name);
% logger.setFilename(logger_file_name)
% logger.setLogLevel(logger.ALL)
% logger.setCommandWindowLevel(logger.ALL)

