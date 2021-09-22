% Oskar Berntsson, 2021
% This script shows an example on how the data analysis could be set up.
%%
close all
clear
scrsz = get(0,'screensize'); % get info on the screensize - might be useful for making plots the "right" size.
%%
% Where to find the data
dataFolder.SAXS = './sampleData/SAXSdata'; % Wgere to find SAXS data
dataFolder.WAXS =  './sampleData/WAXSdata'; % Where to find WAXS data

% Info on the setup
setupInfo.WAXScalib = 'waxscalib.txt'; % Mythen calibration file
setupInfo.detector_readoutrate = 500; % Hz - to get the time-vector correct.
setupInfo.nStepsPerCycle = 2500; % the number of steps per cycle.
setupInfo.dq = 0.005; % the q-step. For rebinning the q-vector.
setupInfo.scaleWAXS = 2.86E-4; % approximately puts SAXS and WAXS data on the same scale

%%
scanNumbers = [642:645]; % Which scans (images) to load

timeShift = 1; % The detector-laser offset, s
imPerCycle = setupInfo.nStepsPerCycle; % the number of images per cycle
replot = false; % just replot the data, not reloading
timeSlices = [-0.5, 0, 0.1, 0.2, 0.3 0.4, 0.5, 1, 2, 3, 4, 5]; % how to divide the time
qSlices = [0.1 0.15];% different q-ranges to look at, for kinetic tracing

plot_dI = true; % whether to plot I or dI
qPower = 0; % plot data as q^qPower*dS
qRange_SAXS = [0.02,0.5]; % show the data in this q-range (1/Å)
qRange_WAXS = [1.4 2.4];
smoothspan = 1; % smoothing
qScale = 'linear'; %'linear' or 'log'
offset = 0.1; % The curves are offset with this value

% SVD
nComp_SVD =1; % number of SVD components to plot
qRange_SVD = [0.03 0.5]; % over which q-range do you want to perform the SVD
qPower_SVD = 1; % what qPower would you like to use for your SVD analysis

monRange = [0.01 0.1]; % over which q-range to monitor drifts in data.

% Outliers and normalization
outlierRange = [2.02 2.12]; % look for outliers in this q-range
outlierLevel = '0.2prcnt'; % 'Xprcnt', 'Xsigma' or 'Xunits'
normRange = [1.45, 1.55]; % normalize to the scattering in this range

% Heat subtraction
do_subtractHeat = true; % false = do not subtract heat, true = do subtract heat
heatData = 'buffer_heating.mat'; % the file with the heat data. The heat data is expected to be stored as a variable named dI_heat with size nq x nHeatComponents
qRange_heat = [1.45 2.4]; % q-range in which to scale heat data to regular data

%% Load the data

% Initialize some variables
I_stack = [];
I_noLaser = [];
I_Laser = [];
dI_stack = [];
dI_all = [];
t_all = [];
scanTag = [];
sumScat = [];

if ~replot
    % Load all the data specified by scannumbers
    for i = 1:length(scanNumbers)
        
        [q, I_, dI_, t_,sumScat_] = rrLoad(scanNumbers(i),normRange,timeShift,outlierLevel,outlierRange, monRange,dataFolder,setupInfo);
        
        dI_stack = [dI_stack, mean(dI_,2)];
        dI_all = [dI_all, dI_];
        
        I_stack = [I_stack, mean(I_.noLaser,2)];
        I_noLaser = [I_noLaser, I_.noLaser];
        I_Laser = [I_Laser, I_.Laser];
        
        sumScat = [sumScat, sumScat_];
        t_all = [t_all,t_];
        scanTag = [scanTag, repmat(scanNumbers(i),1,numel(sumScat_))];
        
    end
    
    % Save data for fast replot
    save('dataTmp.mat','q','I_stack','I_noLaser','I_Laser','dI_all','dI_stack','t_all','scanTag','scanNumbers','sumScat', 'setupInfo','dataFolder')
else
    % Reload previously loaded data, useful for fast replotting
    load('dataTmp.mat')
end

%% Diagnostics figures
figure
set(gcf,'position',[1931         858        1887         250])
hold on
for i = 1:numel(scanNumbers)
    plot(find(scanTag==scanNumbers(i)), sumScat(scanTag==scanNumbers(i)));
end
ah = gca;
ah.XLim = [1 numel(sumScat)];
ah.XTick = 1:setupInfo.nStepsPerCycle:numel(sumScat);
ah.XTickLabel = [];

%xlabel('Image')
ylabel('\SigmaI(a<q<b)')
title('Signal over time')

%% Calculate the average scattering per timepoint
t = unique(t_all);
for i = 1:numel(t)
    dI_avg(:,i) = mean(dI_all(:,t_all==t(i)),2);
    I_avg(:,i) = mean(I_Laser(:,t_all==t(i)),2);
end

%% Subtract heat
if do_subtractHeat
    disp('Will subtract buffer heating')
    load(heatData)
    [dI_avg,ampl_heat] = subtractHeat(q, dI_avg, dI_heat, qRange_heat);
    [dI_stack,~] = subtractHeat(q, dI_stack, dI_heat, qRange_heat);
end

%% Plot the data
% These plots shows some types of analysis that one might want to do
%close all
%offset = 0.02;
if qPower == 0
    %offset = 0.01;
elseif qPower == 1
    %offset = 0.001;
    offset = offset/10;
elseif qPower == 2
    offset = offset/100;
end

fh = figure;
set(fh,'position',[1921 1 1920 1116])


% Regular 1D I vs q type of plot
for i = 1:numel(timeSlices)-1
    dI_tmp(:,i) = mean(dI_avg(:,t>=timeSlices(i)&t<=timeSlices(i+1)),2);
    
    %I_tmp(:,i) = mean(Ilaser(:,t_all>=timeSlices(i)&t_all<=timeSlices(i+1)),2);
    legs{i} = sprintf('%.2f<t<%.2f s',timeSlices(i),timeSlices(i+1));
end

smoothSpan = [];
if plot_dI
    % Plot the SAXS data in one subplot
    subplot(2,4,1)
    stackPlot(q,dI_tmp,qPower,qRange_SAXS,smoothSpan,offset,qScale,legs,'legbox')
    
    % Plot the WAXS data in one subplot
    subplot(2,4,2)
    stackPlot(q,dI_tmp,qPower,qRange_WAXS,smoothSpan,offset,qScale,legs,'legbox')
end
title('Difference scattering at different times')

% Contour plot
subplot(2,4,3)
smoothSpan = [];
if plot_dI
    qt2DPlot(q,t,dI_avg,qPower,qRange_SAXS,smoothSpan)
else
    qt2DPlot(q,t,I_avg,qPower,qRange,smoothSpan)
end
title('2D plot')

% Kinetic traces
subplot(2,4,4)
traceWhat = 'avg';
doNormalize = false;
if plot_dI
    kintracePlot(dI_avg,q,t,traceWhat,qSlices,doNormalize)
else
    kintracePlot(I_avg,q,t,traceWhat,qSlices,doNormalize)
end
title('Kinetic trace')

% Stacked plot
smoothSpan = [];
if plot_dI
    % Plot the SAXS data in one subplot
    subplot(2,4,5)
    stackPlot(q,dI_stack,qPower,qRange_SAXS,smoothSpan,0,qScale,scanNumbers,'legbox')
    
    % Plot the WAXS data in one subplot
    subplot(2,4,6)
    stackPlot(q,dI_stack,qPower,qRange_WAXS,smoothSpan,0,qScale,scanNumbers,'legbox')
else
    % Plot the absolute scattering
    subplot(2,4,5)
    stackPlot(q,I_stack,0,qRange_SAXS,smoothSpan,0,'log',scanNumbers,'legbox')
    set(gca,'yscale','log')
    subplot(2,4,6)
    stackPlot(q,I_stack,0,qRange_WAXS,smoothSpan,0,'linear',scanNumbers,'legbox')
end
title('Stacked plot')

% Plot the absolute scattering
% subplot(2,4,7)
% stackPlot(q,I_stack,0,qRange_SAXS,smoothSpan,0,'log',scanNumbers,'legbox')
% set(gca,'yscale','log')
% subplot(2,4,8)
% stackPlot(q,I_stack,0,qRange_WAXS,smoothSpan,0,'linear',scanNumbers,'legbox')

%
%Peak positions
subplot(2,4,[7 8])
if plot_dI
    peakPosPlot(dI_avg,q,t,[0.02 0.1],qPower,'truff')
else
    peakPosPlot(I_avg,q,t,[0.05 0.11; 0.1 0.2],qPower,'peak')
end
title('Peak or truff position')

%% SVD analysis
if plot_dI
    [U_heat,V_heat] = trxssSVD(q,dI_avg,t,nComp_SVD,qPower_SVD,qRange_heat);
    [U,V] = trxssSVD(q,dI_avg,t,nComp_SVD,qPower_SVD,qRange_SVD);
else
    [U,V] = trxssSVD(q,I_avg,t,nComp_SVD,qPower_SVD,qRange_SVD);
end
fh = gcf;
fh.Position = [76.9673   12.2238   23.9977    7.9904];



