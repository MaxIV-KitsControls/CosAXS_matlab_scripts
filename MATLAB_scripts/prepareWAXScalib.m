function prepareWAXScalib(dataFolder, scanNumber)
% Oskar Berntsson, 2021
%% prepareWAXScalib
% This will help you prepare a calibration file that, via MAIN_PLOT_RR is
% passed to rrLoad. This is required to make the indices form the Mythen
% detector into q-values.
%
% q-values for LAB6 peaks in the WAXS range: 1.512, 2.138
% q-values for Lupolen in the WAXS range: 1.53, 1.69, 2.12

%dataFolder =  '/data/visitors/cosaxs/20200104/2021063008/raw';
close all
if ~isempty(scanNumber)
    WAXSfile = [dataFolder filesep sprintf('mythen_scan_%u_data.hdf5',scanNumber)];
else
    WAXSfile = dataFolder;
end
tmpInfo = h5info(WAXSfile);
nq = tmpInfo.Datasets(1).Dataspace.Size(1);
inds = 1:nq;
I = h5read(WAXSfile,'/data');
I = flipud(I);
I = mean(I,2);

plot(inds,I)
xlabel('Index')
ylabel('Intensity')
%title(sprintf('scan%u',scanNumber))
title(WAXSfile,'interpreter','none')

nPeaks = input('Enter number of peaks: ');
for i = 1:nPeaks
    str = sprintf('Enter index for position of peak %u: ',i);
    peakInds(i,2) = input(str);
    str = sprintf('Enter q-value for peak %u: ',i);
    peakInds(i,1) = input(str);
end

calibFile = input('Enter name of calibration file: ','s');
calibFile = [calibFile '.txt'];
whatToDo = input('Overwrite (o) or append (a) to previous file: ','s');

if ~strcmpi(whatToDo,'o') && ~strcmpi(whatToDo,'a')
    disp('Make a valid choice')
    whatToDo = input('Overwrite (o) or append (a) to previous file: ','s');
end

if strcmpi(whatToDo,'o')
    %writematrix(peakInds,calibFile)
    dlmwrite(calibFile,peakInds)
    
elseif strcmpi(whatToDo,'a')
    tmp = dlmread(calibFile);
    %tmp = readmatrix(calibFile);
    peakInds = [peakInds; tmp];
    peakInds = sortrows(peakInds);
    %writematrix(peakInds,calibFile)
    dlmwrite(calibFile,peakInds)
    
end

p = polyfit(peakInds(:,2),peakInds(:,1),1);
plot(peakInds(:,2),peakInds(:,1),'o',peakInds(:,2),polyval(p,peakInds(:,2)),'r')
xlabel('Index')
ylabel('q-value')
legend('Entered datapoints','Linear fit')

end
