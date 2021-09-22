function [U,V] = trxssSVD(q,dI,t,nComp,qPower,qRange)
%Oskar Berntsson, 2021
%Performs an SVD analysis and plots the results. Returns the nComp first 
%compontents of the U and V vector. Note that U will be weighted by qPower.

fh=figure;
fh.Units = 'centimeters';
fh.Position = [18 15 24 8];

qi = q(q>=qRange(1) & q<=qRange(2));
inData = repmat(q.^qPower,1,size(dI,2)).*dI(:,:);
inData = inData(q>=qRange(1) & q<=qRange(2),:);

[U, S, V] = svd(inData);

subplot(1,3,1)
plot(qi,U(:,1:nComp))
ah = gca;
ah.XLim = qRange;
%ah.XScale = qScale;
title('Left singular vectors')
xlabel('q (Å^{-1})')
ylabel('LSV')

subplot(1,3,2)
plot(t,V(:,1:nComp),'-o')
title('Right singular vectors')
xlabel('Time (s)')
%xlabel('\DeltaTemp. (deg. C)')
%xlabel('T_{final} (deg. C)')
ylabel('Amplitude (arb.)')

subplot(1,3,3)
plot(diag(S)./sum(diag(S)),'-ok')
ah = gca;
ah.XLim = [1, 10];
title('Singular values')
ylabel('Amplitude (arb.)')

U = U(:,1:nComp);
V = V(:,1:nComp);

end

