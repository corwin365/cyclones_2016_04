function LocalDataDir = LocalDataDir()

%identifies which machine we're using, and sets the path to the data store directory accordingly


[~,SystemName] = system('hostname');

%internal system setup details removed for upload to Github. Just set the output to the appropriate path on your system.

return
end