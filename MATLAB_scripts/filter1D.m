function filterOut = filter1D( filterIn, vec1, outlierType,  fignr )
% Oskar Berntsson, 2021
% filterOut = filter1D( filterIn, vec1, outlierType,  fignr )
% applies a histogram filter based on the deviation from the median value
% of any 1D vector (vec1). All which deviate by more than a cutoff value
% are discarded. The outlierType argument can be (with X any number):
%   * 'Xsigma': cutoff is X standard deviations
%   * 'Xprcnt': cutoff is X percent of the median value
%   * 'Xunits': cutoff is absolute X
%
% Alternatively, an explicit range can be specified. For example like
% outlierType = '[1,2]range'

if ~strcmp(outlierType,'none')
    mu = nanmedian(vec1(filterIn));
    
    if strcmp(outlierType(length(outlierType)-4:length(outlierType)),'sigma')
        cutoff = str2num(outlierType(1:length(outlierType)-5))*nanstd(vec1);
        
    elseif strcmp(outlierType(length(outlierType)-4:length(outlierType)),'prcnt')
        cutoff = abs(str2num(outlierType(1:length(outlierType)-5))/100*mu);
        
    elseif strcmp(outlierType(length(outlierType)-4:length(outlierType)),'units')
        cutoff = str2num(outlierType(1:length(outlierType)-5));
        
    elseif strcmp(outlierType(length(outlierType)-4:length(outlierType)),'range')
        range = str2num(outlierType(1:length(outlierType)-5));
        mu = mean(range);
        cutoff = abs(range(2) - range(1)) / 2;
        
    end
    
    bad = find(vec1 < mu - cutoff | vec1 > mu + cutoff);
    
    
    if ~isempty(fignr)
        fh = figure(fignr);
        hist(vec1, length(vec1)/2);
        title(sprintf('outlierType: %s ', outlierType))
        hold on
        YLim = get(gca, 'YLim');
        line([mu - cutoff, mu - cutoff], [YLim(1) YLim(2)],'color','r')
        line([mu + cutoff, mu + cutoff], [YLim(1) YLim(2)],'color','r')
        hold off
        drawnow
        fprintf('    Found %.f outliers of %.f images \n',numel(bad),sum(filterIn));
        if numel(bad) > 0 && numel(bad) < 20
            fprintf('    This was images %s \n',num2str(bad))
        end
    end
    
    filterTemp = true(size(filterIn));
    filterTemp(bad) = false;
    filterOut = filterTemp & filterIn;
else
    filterOut = filterIn;
end


end

