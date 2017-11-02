function fullFilename = saveResults(results,dirname)

if nargin < 2 || isempty(dirname)
    dirname = 'Results';
end

geometry = results(1).geometry;

dimStr = sprintf('%dd',geometry.dim);
geometryStr = geometry.name;

specialSettingStr = geometry.specialSetting;

if strcmpi(specialSettingStr,'none')
    specialSettingStr = [];
end

diagJumpStr = [];
if geometry.diagJumps == 1
    diagJumpStr = 'diagJumps';
elseif geometry.diagJumps == 2
    diagJumpStr = 'diagJumpsCorrected';
end

c = clock;
dateStr = sprintf('%.2d_%.2d_%.2d_%.2d_%.2d_%.2d',c(1),c(2),c(3),c(4),c(5),round(c(6)));

filename = [dimStr '_' geometryStr '_' specialSettingStr '_' diagJumpStr '_' dateStr];
filename = strrep(filename,'__','_');
filename = strrep(filename,'__','_');
filename = strrep(filename,'__','_');

fullFilename = [dirname '/' filename];

if ~exist(dirname,'dir')
    mkdir(dirname);
end

save(fullFilename,'results');
fprintf('Results saved in %s.\n',fullFilename);

end