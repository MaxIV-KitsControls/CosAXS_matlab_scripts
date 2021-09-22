function [] = stackPlot(q,dI,qPower,qRange,smoothSpan,offset,qScale,legs,legPos)
% Oskar Berntsson, 2021
%Makes a regular 1D plot. Apart from the general data visualization this is 
%also useful for monitoring scan to scan variation,concentration series, 
%or power titration. If legs are supplied as a
%numerical vector stackPlot assumes that you want to monitor scan to scan 
%variation, if you want to visualize a power titration, provide appropriate 
%legs as a cell array
x = q;
if isempty(smoothSpan)
    smoothSpan = 1;
end

y = repmat(q.^qPower,1,size(dI,2)).*dI;
for i = 1:size(y,2)
   %y(:,i) = smooth(y(:,i),smoothSpan) + (1-i)*offset;
   y(:,i) = medfilt1(y(:,i),smoothSpan) + (1-i)*offset; 
end

plot(x,y);


ah = gca;
ah.XLim = qRange;
ah.XScale = qScale;

xlabel('q (Å^{-1})')
if qPower == 0
    ylabel('\DeltaI(q,t) (arb.)')
elseif qPower == 1
    ylabel('q\DeltaI(q,t) (arb.)')
else
    ylabel(['q^' num2str(qPower) '\DeltaI(q,t) (arb.)'])
end

if ~isempty(legs)
    if isnumeric(legs)
        for i = 1:numel(legs)
            leg{i} = sprintf('Scan%u',legs(i));
        end
    else
        leg = legs;
    end
    
    if strcmpi(legPos,'legbox')
        legend(leg,'location','best')
    elseif strcmpi(legPos,'online')
        for i = 1:numel(legs)            
            %text(qRange(2)*0.7,0.5*offset+(1-i)*offset,legs{i},'BackgroundColor','white');
            text(0.6,0.5*offset+(1-i)*offset,legs{i});
        end
    end
end

grid on
box on


end

