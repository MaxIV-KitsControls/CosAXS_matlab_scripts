function [] = peakPosPlot(dI,q,t,qRange,qPower,peakOrtruff)
% Oskar Berntsson, 2021
%finds the position of a peak or truff and plots the position in q vs time.


dI = repmat(q.^qPower,1,size(dI,2)).*dI;
%dI = abs(dI);
y = medfilt1(dI,10);

for i = 1:size(qRange,1)
[dI_tmp, ~] = qCut(dI,q,qRange(i,:));
[y_tmp, q_tmp] = qCut(y,q,qRange(i,:));
leg{i} = sprintf('%.2f<q<%.2f Å^{-1}',qRange(i,1),qRange(i,2));
for ii = 1:size(dI_tmp,2)
        
%            plot(q_tmp,dI_tmp(:,ii),'.k',q_tmp,y_tmp(:,ii),'r')
%            drawnow
%            pause(0.01)       
        
        if strcmpi(peakOrtruff,'peak')
            inds = q_tmp(y_tmp(:,ii)==max(y_tmp(:,ii)));
            peakPos(i,ii) = inds(1);
            ylab = 'Peak position (Å^{-1})';
        elseif strcmpi(peakOrtruff,'truff')
            inds = q_tmp(y_tmp(:,ii)==min(y_tmp(:,ii)));
            peakPos(i,ii) = inds(1);
            ylab = 'Truff position (Å^{-1})';
        end
       if t(ii) <0
           peakPos(i,ii) = NaN;
       end
end


plot(t,peakPos,'.')
xlim([min(t) max(t)])
xlabel('Time (s)')
ylabel(ylab)
legend(leg,'location','best')

grid on
box on


end

