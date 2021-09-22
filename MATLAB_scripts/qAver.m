function Iout = qAver( I, q, qRange )
% Oskar Berntsson, 2021
% Iout = qAver( I, q, qRange )
% averages the scattering data in I over the specified q range, and returns
% a row vector holding the average for each image.

%Iout = nanmean(I(q >= min(qRange) & q <= max(qRange),:),1);
Iout = mean(I(q >= min(qRange) & q <= max(qRange),:),1,'omitnan');

end
