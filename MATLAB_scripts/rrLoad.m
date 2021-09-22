function [q, I, dI, t,sumScat]=rrLoad(scanNumber,normRange,timeShift,outlierLevel,outlierRange, monRange, dataFolder,setupInfo)
% Oskar Berntsson, 2021
%% INPUT
% scanNumber:    The scan number
% normRange:     The q-range over which to normalize
% timeShift:     The delay time between lasers and X-rays. Used to calculate timepoints
% outlierLevel:  By which measure are outliers judged
% outlierTange:  The q-range over which to look for outliers
% dataFolder:    structure with fields SAXS and WAXS which tells where the data is
% setupInfo:     structure with fields about the setup (wavelength, detector distance, etc.)

%% OUTPUT
% q
% I
% dI
% t
% sumScat

%% SUBFUNCTIONS
% qSum
% qAver
% cycleReject
%% SETTINGS
nStepsPerCycle = setupInfo.nStepsPerCycle;
t = (((1:nStepsPerCycle)-1)*1/setupInfo.detector_readoutrate)-timeShift; % set up the time vector

% Initialize some variables
whichStep = 1;
nGood = 0;
nGoodON = 0;
nGoodOFF = 0;
dI = 0;
sumScat = [];
I_.Laser = 0;
I_.noLaser = 0;
I_all = [];
%% LOAD DATA
fprintf('SCAN %u\n',scanNumber) % print which data is to be loaded
%fprintf('    Loading data\n') % print which data is to be loaded

% Define and get some info from the datafiles
WAXSfile = [dataFolder.WAXS filesep sprintf('mythen_scan_%u_data.hdf5',scanNumber)];
SAXSfile = [dataFolder.SAXS filesep sprintf('eiger_scan_%u_data.hdf5',scanNumber)];

% Get some info from the files
tmpInfo = h5info(WAXSfile);
nqWAXS = tmpInfo.Datasets(1).Dataspace.Size(1);
nCycles(1) = floor(tmpInfo.Datasets(1).Dataspace.Size(2)/nStepsPerCycle); % this means that we only account for complete cycles
nData(1) = tmpInfo.Datasets(1).Dataspace.Size(2);

tmpInfo = h5info(SAXSfile);
nqSAXS = tmpInfo.Datasets(1).Dataspace.Size(1);
nCycles(2) = floor(tmpInfo.Datasets(1).Dataspace.Size(2)/nStepsPerCycle);
nData(2) = tmpInfo.Datasets(1).Dataspace.Size(2);

% Check if number of frames in SAXS and WAXS files are the same
if nData(1) ~= nData(2)
    droppedFrames = true;
    fprintf('Found %u SAXS curves and %u WAXS curves\n',nData(2),nData(1));
    fprintf('Will add NaN for missing curves and disregard these for the averages\n')
    
    % Sometimes the filewriter drop frames from the Eiger. These are only a few frames (about 0.02%).
    % They will likely not affect the data unless they are "unfortunately" placed.
    % This part of the code figures out what data has been dropped so that
    % this can be replaced by NaNs later on //OB Aug.2021
    
    frameInd = h5read([dataFolder.RAW filesep sprintf('eiger_scan_%u_data.hdf5',scanNumber)],'/entry/data/meta'); % load the frame count
    frameInd = single(frameInd(1,:))+1; % adding +1 to agree with matlab indexing
    frameInd_ = diff(frameInd); % check which frames are not consequtive
    frameInd_ = [1 frameInd_]; % make indices agree with frameInd
    dropInd = find(frameInd_>1); % where were frames dropped, starting position, matlab indices
    dropLen = frameInd_(dropInd)-1; % how many frames where dropped  
    
    % This loop adds NaN where there where dropped frames. 
    for i = 1:numel(dropInd)
        tmp = [frameInd(1:dropInd(i)-1), NaN(1,dropLen(i)), frameInd(dropInd(i):end)];
        if i<numel(dropInd)
            dropInd(i+1) = dropInd(i+1)+dropLen(i);
        end
        frameInd = tmp;
    end
    
    missingFrame = find(isnan(frameInd)); % which frames are missing    
    dropCycle = ceil(missingFrame/nStepsPerCycle); % find in which cycle there are dropped frames
    dropInd = mod(missingFrame,nStepsPerCycle); % find where within a cycle the drop is    
    dropCycleUnique = unique(dropCycle);
    
    for i = 1:numel(dropCycleUnique)
        droppedFramesPerCycle(i) = sum(dropCycle == dropCycleUnique(i));
    end
    
    droppedFramesTotal = 0;
    
else
    droppedFrames = false;
end
nCycles = nCycles(1); % The Mythen data is always correct, i.e. no frames are dropped.

WAXScount = [nqWAXS nStepsPerCycle]; % how much should be read from the WAXS file, per cycle
SAXScount = [nqSAXS nStepsPerCycle]; % how much should be read from the SAXS file, per cycle

qSAXS = h5read(SAXSfile,'/q'); % read the q-vector from SAXS file
qSAXS = qSAXS/10; % from 1/nm to 1/A

% Calibrate WAXS range
%peak_inds = [439 449 587 949 960]; % Lupolen+LAB6, Scans 438 and 442
%q_vals_standard = [1.512 1.53 1.69 2.12 2.138]; % Lupolen+LAB6,
%tmp = readmatrix(setupInfo.WAXScalib); % load a calibration file
tmp = dlmread(setupInfo.WAXScalib); % load a calibration file
q_vals_standard = tmp(:,1);
peak_inds = tmp(:,2);
p = polyfit(peak_inds,q_vals_standard,1);
r_inds = 1:nqWAXS;
qWAXS = polyval(p,r_inds)';

% Setup the q-vector
qlims = [min(qSAXS) max(qWAXS)];
q = [qlims(1):setupInfo.dq:qlims(2)]';

% To make the rebinning work, if not rounded - it sometimes thinks
% qtmp~=q(i) because of precision "floatsequal"
qi = round(q./setupInfo.dq);
qtmp = round([qSAXS; qWAXS]./setupInfo.dq);

fprintf('    Found %u cycles\n',nCycles);
if ~outlierLevel | strcmpi('none',outlierLevel) | isempty(outlierRange)
    disp('    No outlier rejection!')
else
    disp('    Will perform outlier rejection!')
end

if ~isempty(normRange)
    fprintf('    Will normalize data at %.2f<q<%.2f\n',normRange(1),normRange(2));
else
    disp('    No normalization')
end

for cycleNumber = 1:nCycles
    %fprintf('    Cycle %u/%u\n',cycleNumber,nCycles);
    
    % Load WAXS data
    WAXSstart = [1 (cycleNumber-1)*nStepsPerCycle+1];
    IWAXS = h5read(WAXSfile,'/data',WAXSstart,WAXScount);
    IWAXS = flipud(IWAXS); % because Mythen is mounted upside down
    
    % Load SAXS data
    if droppedFrames
        % If there are dropped frames, this part of the code loads the
        % appropriate data for each cycle and adds NaN for missing data. 
        % The NaNs will be disregarded/not counted in the averaging process later on.
        % //OB Aug. 2021        
        droppedFramesThisCycle = droppedFramesPerCycle(dropCycleUnique==cycleNumber);
        if isempty(droppedFramesThisCycle)
            droppedFramesThisCycle = 0;
        end
        
        % Adjust the from and to so that only data for this cycle are read
        SAXSstart = double([1 (cycleNumber-1)*nStepsPerCycle+1-droppedFramesTotal]);
        droppedFramesTotal = droppedFramesTotal + droppedFramesThisCycle;         
        SAXScount = double([nqSAXS nStepsPerCycle - droppedFramesThisCycle]); % how much should be read from the SAXS file, per cycle
        
        ISAXS = h5read(SAXSfile,'/I',SAXSstart,SAXScount);
        
        if droppedFramesThisCycle>0
            
            dropInd_tmp = dropInd(dropCycle == cycleNumber);
            % This loop adds NaN where there were dropped frames
            for i = 1:numel(dropInd_tmp)
                tmp = [ISAXS(:,1:dropInd_tmp(i)-1), NaN(nqSAXS,1), ISAXS(:,dropInd_tmp(i):end)];
                if i<numel(dropInd_tmp)
                    dropInd_tmp(i+1) = dropInd_tmp(i+1)+1;
                end
                ISAXS = tmp;
            end
            
        end       
        
    else
        % If no frames are dropped - proceed as usual.
        SAXSstart = [1 (cycleNumber-1)*nStepsPerCycle+1];
        ISAXS = h5read(SAXSfile,'/I',SAXSstart,SAXScount);
    end
    
    % Rebin data
    Itmp = [ISAXS;setupInfo.scaleWAXS*single(IWAXS)];
    I = zeros(numel(q),nStepsPerCycle);
    for i = 1:numel(qi)
        I(i,:) = mean(Itmp(qtmp==qi(i),:),1);
    end
    
    % Cut range that cannot be trusted for whatever reason
    I(q>0.5&q<1.4,:) = NaN;
    
    % Sum the scattering in some q-region(s). This can be useful for
    % diagnostic purposes.
    for i = 1:size(monRange,1)
        sumScat = cat(2,sumScat,qSum(I,q,monRange(i,:)));
    end
    
    %fprintf('    Data loaded\n'); % update on progress
    
    %% OUTLIER REJECTION
    %fprintf('    Outlier rejection\n'); % update on progress
    if ~outlierLevel | strcmpi('none',outlierLevel) | isempty(outlierRange)
        disp('No outlier rejection!')
    else
        vec1 = qAver(I,q,outlierRange);
        % Laser heating would result in a major change of the scattering
        % signal - and might be detected as an outlier cycle.
        % But if you monitor the region around the isosbestic
        % points at 1.5 or 2.07 i.e. outlierRange = [1.45 1.55] or [2.02 2.12]
        % this should work
        h = cycleReject(vec1,outlierLevel); % h=0 means data is "normally" distributed, h=1 means data is not
        
        if h == 1
            isBadCycle = true;
            fprintf('Possibly bad data in scan%u, cycle %u\n',scanNumber, cycleNumber)
            histogram(vec1);
            ah = gca;
            line([mean(vec1),mean(vec1)],[ah.YLim(1) ah.YLim(2)],'Color','red','LineStyle','--');
            line([median(vec1),median(vec1)],[ah.YLim(1) ah.YLim(2)],'Color','blue','LineStyle','--');
            title(sprintf('Step %u in scan%u',cycleNumber,scanNumber))
            
            pause(1)
        else
            isBadCycle = false;
        end
        
    end
    
    %% NORMALIZATION
    %fprintf('    Normalization\n'); % update on progress
    if ~isempty(normRange)
        I = normalize( q, I, normRange );
    else
        disp('No normalization')
    end
    
    %% Figure out cycles and calculate differences
    % odd cycles are always laser off
    
    if whichStep == 1
        % First cycle should be laser off
        Ibefore = I;
        whichStep = 2;
        
        if ~isBadCycle
            nGoodOFF = nGoodOFF+1;
            
            % OLD
            %I_.noLaser = I_.noLaser+Ibefore; % this is to store all non-subtracted data from "laser off" cycles
            
            % New: this is to accomodate that frames dropped by the
            % filewriter are replaced by NaNs //OB Aug.2021
            if nGoodOFF == 1
                I_.noLaser = Ibefore;
            else
                I_.noLaser = 2*mean(cat(3,(nGoodOFF-1)/nGoodOFF*I_.noLaser,1/nGoodOFF*Ibefore),3,'omitnan'); % this is to store all non-subtracted data from "laser off" cycles
            end
            
        end
        
    elseif whichStep == 2
        % even cycles should be laser on
        Ilaser = I;
        whichStep = 3;
        
        if ~isBadCycle
            nGoodON = nGoodON+1;
            
            % OLD
            %I_.Laser = I_.Laser+I; % this is to store all non-subtracted data from "laser on" cycles
            
            % New: this is to accomodate that frames dropped by the
            % filewriter are replaced by NaNs //OB Aug.2021
            if nGoodON == 1
                I_.Laser = Ibefore;
            else
                I_.Laser = 2*mean(cat(3,(nGoodON-1)/nGoodON*I_.Laser,1/nGoodON*Ilaser),3,'omitnan'); % this is to store all non-subtracted data from "laser on" cycles
            end
            
        end
        
    elseif whichStep == 3
        Iafter = I;
        
        % OLD
        %dI_tmp = Ilaser - (Ibefore+Iafter)./2;
        
        % New: this is to accomodate that frames dropped by the
        % filewriter are replaced by NaNs //OB Aug.2021
        dI_tmp = Ilaser - mean(cat(3,Ibefore,Iafter),3,'omitnan');
        
        if ~isBadCycle
            nGood = nGood+1;
            % OLD
            %dI = dI + dI_tmp;
            
            % New: this is to accomodate that frames dropped by the
            % filewriter are replaced by NaNs //OB Aug.2021
            if nGood == 1
                dI = dI_tmp;
            else
                dI = 2*mean(cat(3,(nGood-1)/nGood*dI,1/nGood*dI_tmp),3,'omitnan');
            end
            
            
        end
        
        if ~isBadCycle
            nGoodOFF = nGoodOFF+1;
            % OLD
            %I_.noLaser = I_.noLaser+Iafter; % this is to store all non-subtracted data from "laser off" cycles
            
            % New: this is to accomodate that frames dropped by the
            % filewriter are replaced by NaNs //OB Aug.2021
            if nGoodOFF == 1
                I_.noLaser = Iafter;
            else
                I_.noLaser = 2*mean(cat(3,(nGoodOFF-1)/nGoodOFF*I_.noLaser,1/nGoodOFF*Iafter),3,'omitnan'); % this is to store all non-subtracted data from "laser off" cycles
            end
            
        end
        
        Ibefore = I;
        whichStep = 2;
        
    end
    
    % These lines would concatenate all cycles. Might be useful for some
    % diagnostics, but will quickly become massive amounts of data, and
    % slow to process
    %I_all = [I_all, I];
    %dI_all = dI_all, dI]
    
end

clear I
fprintf('SCAN %u, %u Good difference scattering datasets\n',scanNumber,nGood) % print which data was loaded

% As of changes made Aug. 2021 the averageing is already done at an earlier
% stage.
%I.Laser = I_.Laser./nGoodON;
%I.noLaser = I_.noLaser./nGoodOFF;
%dI = dI./nGood;

I.Laser = I_.Laser;
I.noLaser = I_.noLaser;

end

%% SUBFUNCTIONS
% These are functions utilized by rrLoad. They may also exist as
% independent functions that can be used at other places in analysis, but
% rrLoad will use the versions that are listed here.

function Iout = qAver( I, q, qRange )
% Iout = qAver( I, q, qRange )
% averages the scattering data in I over the specified q range, and returns
% a row vector holding the average for each image.

%Iout = nanmean(I(q >= min(qRange) & q <= max(qRange),:),1);
Iout = mean(I(q >= min(qRange) & q <= max(qRange),:),1,'omitnan');

end

function Iout = qSum( I, q, qRange )
% Iout = qSum( I, q, qRange )
% identical to qAver(), but returns a sum instead.

%Iout = nansum(I(q >= min(qRange) & q <= max(qRange),:),1);
Iout = sum(I(q >= min(qRange) & q <= max(qRange),:),1,'omitnan');

end

function [h] = cycleReject(vec1, outlierLevel)
% Compares the mean and median of vec1 and considers the cycle to be an
% outlier if the two are more different than what is given by outlierLevel.
%%
if strcmp(outlierLevel(length(outlierLevel)-4:length(outlierLevel)),'sigma')
    datComp = mean(vec1)-median(vec1);
    if abs(datComp) > str2num(outlierLevel(1:length(outlierLevel)-5))*nanstd(vec1)
        h = 1;
    else
        h = 0;
    end
    
elseif strcmp(outlierLevel(length(outlierLevel)-4:length(outlierLevel)),'prcnt')
    datFrac = mean(vec1)/median(vec1);
    if abs(1-datFrac) > str2num(outlierLevel(1:length(outlierLevel)-5))/100
        h = 1;
    else
        h = 0;
    end
    
elseif strcmp(outlierLevel(length(outlierLevel)-4:length(outlierLevel)),'units')
    datComp = mean(vec1)-median(vec1);
    
    if abs(datComp) > str2num(outlierLevel(1:length(outlierLevel)-5))
        h = 1;
    else
        h = 0;
    end
end
end
