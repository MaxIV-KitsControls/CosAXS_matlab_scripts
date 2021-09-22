% Oskar Berntsson, 2021
% A simple load data and plot script. This doesn't perform any outlier
% rejection, normalization or calculation of difference scattering.
% It simply loads and avererages data on a per scan basis.
close all
clear
%%
dataFolder.SAXS = './sampleData/SAXSdata'; % SAXS data are found at...
dataFolder.WAXS = './sampleData/WAXSdata'; % WAXS data are found at...
WAXScalib = 'waxscalib.txt'; % Mythen calibration file

dq=0.005; % The desired q-spacing
loadSAXS = true; % Load SAXS data (or not)
loadWAXS = true; % Load WAXS data (or not)
scanNumbers = [519 520]; % which scans to load
      
%% Load the data

% Calibrate WAXS range
%tmp = readmatrix(WAXScalib); % load a calibration file
tmp = dlmread(WAXScalib); % load a calibration file
q_vals_standard = tmp(:,1);
peak_inds = tmp(:,2);
p = polyfit(peak_inds,q_vals_standard,1);


IWAXS = [];
ISAXS = [];
for i = 1:numel(scanNumbers)
    scannumber = scanNumbers(i);
    legs{i} = sprintf('scan%u',scannumber);
    
    if loadWAXS
        WAXSfile = [dataFolder.WAXS filesep sprintf('mythen_scan_%u_data.hdf5',scannumber)];
        tmpInfo = h5info(WAXSfile);
        nqWAXS = tmpInfo.Datasets(1).Dataspace.Size(1);        
        
        if i == 1
            r_inds = 1:nqWAXS;
            qWAXS = polyval(p,r_inds)';            
            
            % For rebining q
            qWAXS = dq*round(qWAXS/dq);
            qi = unique(qWAXS);
            q_tmp = qWAXS;
            qWAXS = qi;
        end
                
        
        % Load WAXS data
        I_ = h5read(WAXSfile,'/data');
        I_ = mean(I_,2);
        I_ = flipud(I_);
        
        % Rebin q-range for WAXS
        I_tmp = zeros(numel(qi),size(I_,2));
        for i = 1:numel(qi)
            I_tmp(i,:) = mean(I_(q_tmp==qi(i),:));
        end
        IWAXS = [IWAXS, I_tmp];
        
    end
    
    if loadSAXS
        SAXSfile = [dataFolder.SAXS filesep sprintf('eiger_scan_%u_data.hdf5',scannumber)];
        tmpInfo = h5info(SAXSfile);
        nqSAXS = tmpInfo.Datasets(1).Dataspace.Size(1);
        
        qSAXS = h5read(SAXSfile,'/q');
        qSAXS = qSAXS/10; % from 1/nm to 1/A
        
        % Load SAXS data
        I_tmp = h5read(SAXSfile,'/I');
        intensities = qSum(I_tmp,qSAXS,[0.02 0.5]);
        I_tmp = mean(I_tmp,2);
        
        ISAXS = [ISAXS, I_tmp];
    end
    
    
    
end


%% Plot the data

fh = figure;
fh.Position = [1447         693        1109         406];

if loadSAXS
    subplot(1,2,1)
    plot(qSAXS,ISAXS)
    ah = gca;
    ah.YScale = 'log';
    ah.XLim = [0 max(qSAXS)];    
    xlabel('q [1/Å]')
    ylabel('I(q) [arb]')
    
    legend(legs)
    title('SAXS')
    
end

if loadWAXS
    subplot(1,2,2)
    plot(qWAXS,IWAXS)
    ah = gca;
    %ah.YScale = 'log';
    ah.XLim = [min(qWAXS) max(qWAXS)];
    xlabel('q [1/Å]')
    ylabel('I(q) [arb]')   

    legend(legs)
    title('WAXS')
end


