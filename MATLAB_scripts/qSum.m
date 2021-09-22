function Iout = qSum( I, q, qRange )
% Oskar Berntsson, 2021
% Iout = qSum( I, q, qRange )
% identical to qAver(), but returns a sum instead.

%Iout = nansum(I(q >= min(qRange) & q <= max(qRange),:),1);
Iout = sum(I(q >= min(qRange) & q <= max(qRange),:),1,'omitnan');

end

