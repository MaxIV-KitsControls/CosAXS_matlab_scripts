function [] = kintracePlot(dI,q,t,traceWhat,qSlices,doNormalize)
% Oskar Berntsson, 2021
%Makes a plot showing the change in intensity (sum or average) over a
%certain q-region over time.

for i = 1:size(qSlices,1)
    if strcmpi(traceWhat,'avg')
        kinTrace(i,:) = qAver(dI,q,qSlices(i,:));
    elseif strcmpi(traceWhat,'sum')
        kinTrace(i,:) = qSum(dI,q,qSlices(i,:));
    end
    leg{i} = sprintf('%.2f<q<%.2f Å^{-1}',qSlices(i,1),qSlices(i,2));
end

if doNormalize
    kinTrace = kinTrace ./ repmat(nanmean(kinTrace,2),1,size(kinTrace,2));    
end


plot(t,kinTrace)
xlim([min(t) max(t)])
xlabel('Time (s)')
ylabel('Signal (arb.)')
if doNormalize
    ylabel('Normalized signal')
end
legend(leg,'location','best')
grid on
box on

end

