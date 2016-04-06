clear all

%trim data to make it more tractable for fixed-height analyses
%Corwin Wright, 04/APR/2016

MaxTime = 15;
z = 6;

SABER  = load([LocalDataDir,'/corwin/storms/results/storms_SABER.mat']);

Good = find(SABER.Results.StormIDs(5,:) < MaxTime);
SABER.Results.StormIDs = SABER.Results.StormIDs(:,Good);
SABER.Results.Results  = SABER.Results.Results( z,:,Good);
SABER.Results.Results  = squeeze(SABER.Results.Results);
save([LocalDataDir,'/corwin/storms/results/storms_SABER_atz.mat'],'SABER')


HIRDLS  = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS.mat']);
Good = find(HIRDLS.Results.StormIDs(5,:) < MaxTime);
HIRDLS.Results.StormIDs = HIRDLS.Results.StormIDs(:,Good);
HIRDLS.Results.Results  = HIRDLS.Results.Results( z,:,Good);
HIRDLS.Results.Results  = squeeze(HIRDLS.Results.Results);
save([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_atz.mat'],'HIRDLS')

MLS  = load([LocalDataDir,'/corwin/storms/results/storms_MLS.mat']);
Good = find(MLS.Results.StormIDs(5,:) < MaxTime);
MLS.Results.StormIDs = MLS.Results.StormIDs(:,Good);
MLS.Results.Results  = MLS.Results.Results( z,:,Good);
MLS.Results.Results  = squeeze(MLS.Results.Results);
save([LocalDataDir,'/corwin/storms/results/storms_MLS_atz.mat'],'MLS')