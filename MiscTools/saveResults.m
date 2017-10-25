function fullFilename = saveResults(results)

dirname = sprintf('Results_%dd_%s',results(1).geometry.dim,results(1).geometry.name);

if results(1).geometry.diagJumps == 1
    dirname = [dirname '_diagJumps'];
elseif results(1).geometry.diagJumps == 2
    dirname = [dirname '_diagJumpsCorrected'];
end

filename = sprintf('results_%s',myClock(6));

fullFilename = [dirname '/' filename];

if ~exist(dirname,'dir')
    mkdir(dirname);
end

save(fullFilename,'results');
fprintf('Results saved in %s.\n',fullFilename);

end