clear all

%trim data to make it more tractable for fixed-range analyses
%Corwin Wright, 04/APR/2016

MaxTime = 15; %days
MaxDist = 300e3; %km

SABER  = load([LocalDataDir,'/corwin/storms2/results/storms_SABER.mat']);
Good = find(SABER.Results.StormIDs(6,:) < MaxDist ...
          & SABER.Results.StormIDs(5,:) < MaxTime);
SABER.Results.StormIDs = SABER.Results.StormIDs(:,Good);
SABER.Results.Results  = SABER.Results.Results( :,:,Good);
save([LocalDataDir,'/corwin/storms2/results/storms_SABER_close.mat'],'SABER')

% 
% HIRDLS  = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS.mat']);
% Good = find(HIRDLS.Results.StormIDs(6,:) < MaxDist ...
%           & HIRDLS.Results.StormIDs(5,:) < MaxTime);
% HIRDLS.Results.StormIDs = HIRDLS.Results.StormIDs(:,Good);
% HIRDLS.Results.Results  = HIRDLS.Results.Results( :,:,Good);
% save([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_close.mat'],'HIRDLS')
% 
% MLS  = load([LocalDataDir,'/corwin/storms/results/storms_MLS.mat']);
% Good = find(MLS.Results.StormIDs(6,:) < MaxDist ...
%           & MLS.Results.StormIDs(5,:) < MaxTime);
% MLS.Results.StormIDs = MLS.Results.StormIDs(:,Good);
% MLS.Results.Results  = MLS.Results.Results( :,:,Good);
% save([LocalDataDir,'/corwin/storms/results/storms_MLS_close.mat'],'MLS')