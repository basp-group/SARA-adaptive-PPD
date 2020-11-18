 function output =util_time_based_block_sp_ar(time,param)

if ~isfield(param,'verbose')
    verbose =0;
else
    verbose =param.verbose;
end
if ~isfield(param,'upperTolBlkSz') %defining the limits of accepted block sizes
    upperTolBlkSz =0.5;%%
else
    upperTolBlkSz = param.upperTolBlkSz; %should be in [0 0.5]
end
if ~isfield(param,'lowerTolBlkSz') %defining the limits of accepted block sizes
    lowerTolBlkSz = 0.5 ;%%
else
    lowerTolBlkSz = param.lowerTolBlkSz; %should be in [0 0.5]
end
if ~isfield(param,'tolBlkSzPc') %tolerance on the limits ->might be absorbed in the formers
    tolBlkSzPc = 0.01;%%
else
    tolBlkSzPc = param.tolBlkSzPc; %should be in [0 0.05]
end
output.block = [] ;

%% init
param.pos=[0; param.pos(:)];
%% block size interval
%% block size interval
if param.size<param.pos(end)
   minBlkSz = floor(param.size * (1-lowerTolBlkSz));
   maxBlkSz = floor(param.size * (1+upperTolBlkSz));
else
    minBlkSz =param.size ;
    maxBlkSz =param.size ;
end

fprintf("Block size should be in [%d %d]\n",(1-tolBlkSzPc)*minBlkSz,(1+tolBlkSzPc)*maxBlkSz)

%% partionning
for j =1:length(param.pos)-1
    Blocks = [];
    flag_further_blocking = 1;
    flag_scan_merge = 0;
    dataSlice = time(1+sum(param.pos(1:j)) : sum(param.pos(1:j+1)));
    dataSliceSz =param.pos(j+1) ;
    %  approx. Number & Size of Blocks
    approxBlkNbr =  max(1,floor(dataSliceSz/param.size)) ;
    if floor(dataSliceSz/(approxBlkNbr+1)) >= floor((1-tolBlkSzPc)*minBlkSz)
        approxBlkNbr = approxBlkNbr+1;
    end
    
    BlkSZ = max(floor(dataSliceSz/approxBlkNbr),param.size);
    fprintf("\nINFO: Data slice sz %d",dataSliceSz)
    
    %     fprintf("\nINFO: Data slice of sz %d--> expected nbr and av. sz of blocks %d and %d",dataSliceSz,approxBlkNbr,BlkSZ)
    if approxBlkNbr==1
        Blocks = dataSliceSz;
    else
        %% get snapshot ids and size in the indexed vect. snapshotSzs
        snapshotSzs =[];
        incSnapIdx = 1;
        incSnapSz  = 1;
        for  idBlk = 2 : dataSliceSz
            if dataSlice(idBlk) ~= dataSlice(idBlk-1) %new snapshot
                snapshotSzs(incSnapIdx)= idBlk - incSnapSz;
                incSnapSz =  idBlk;
                incSnapIdx=incSnapIdx+1;
            end
        end
        snapshotSzs(incSnapIdx) = dataSliceSz - incSnapSz+1;
        snapshotsNbr = length(snapshotSzs);
        
        %% identify time scans (by their starting snapshot) if any
        % alternatively read them directly from MS
        diffTimeSnap = zeros(1,snapshotsNbr - 1);
        dummy =0;
        for idBlk=1:snapshotsNbr - 1
            dummy = dummy + snapshotSzs(idBlk);
            currentSnapshotTime = dataSlice(dummy);
            nextSnapshotTime    = dataSlice(dummy+1);
            diffTimeSnap(idBlk)  = nextSnapshotTime -currentSnapshotTime;
        end
        scanIdx=find((diffTimeSnap)>3*median(abs(diffTimeSnap)));
        
        if isempty(scanIdx)
            %% No time scans identified
            %fprintf('\nINFO:No time scans identified')
            %fprintf('\n--> splitting by equal nbr of snapshots')
            approxNbrSnapshotsPerScan =  floor(snapshotsNbr/approxBlkNbr);
            %creating vect. of scans i.e. collection of snapshots (vect. with idx of first
            %snapshot in the scan)
            scanStartingSnap=1:approxNbrSnapshotsPerScan:snapshotsNbr;
            scanStartingSnap(end) = scanStartingSnap(end) + mod(snapshotsNbr,approxBlkNbr);% add remaining snapshots to last scan
            scanEndingSnap = scanStartingSnap(2:end)-1 ;
            scanEndingSnap(end+1)=snapshotsNbr;
            %/ partitionning
            idBlk =1;
            while  idBlk <= approxBlkNbr
                snapLast = scanEndingSnap(idBlk);
                snapFirst= scanStartingSnap(idBlk);
                slice=sum(snapshotSzs(snapFirst:snapLast));
                idBlk=idBlk+1;
                Blocks=[Blocks; slice];
            end
            %handling remaining slice
            if sum(Blocks)~= dataSliceSz
                meanBlkSz=mean(Blocks);
                lastSlice =  dataSliceSz - sum(Blocks);
                if lastSlice >= floor((1-tolBlkSzPc) * min(meanBlkSz,minBlkSz))
                    Blocks=[Blocks(:) ; lastSlice];%adding remainig slice to last block if small
                else
                    Blocks(end)=Blocks(end)+lastSlice ;%creating a new block with remaining slice if big
                end
            end
            %/
            
        else
            %% time scans detected
            %fprintf('\nINFO: %d Time scans detected',length(scanIdx)+1)
            initScanStartingSnap = [1 scanIdx+1];
            %fprintf('--> Initial blocking based on scans')
            %/ init partitionning
            initBlocks =zeros(length(initScanStartingSnap),1);
            for idBlk=1:length(initScanStartingSnap)-1
                initBlocks(idBlk)=sum(snapshotSzs(initScanStartingSnap(idBlk):initScanStartingSnap(idBlk+1)-1));
            end
            initBlocks(length(initScanStartingSnap))=sum(snapshotSzs(initScanStartingSnap(end):end));
            %/
            meanScanSz  = floor(mean(initBlocks));            
            if  meanScanSz < floor(minBlkSz*(1-tolBlkSzPc))
                %/merging scans: snapshots step
                %/ rules are such that  lower bound on the size of the
                %block is a priority
                
                flag_scan_merge =1;
                %fprintf('\nINFO: average scan size is %d --> attempt to merge',meanScanSz)
                nbrScans2Merge = 1+floor(BlkSZ/meanScanSz);
                %/ adjusting snapshots step
                if  meanScanSz*(nbrScans2Merge+1)< maxBlkSz
                    nbrScans2Merge =nbrScans2Merge+1;
                    
                elseif (meanScanSz*nbrScans2Merge  > maxBlkSz*(1+tolBlkSzPc))&& ...
                        (meanScanSz*(nbrScans2Merge-1) >= minBlkSz*(1-tolBlkSzPc) )
                    nbrScans2Merge =nbrScans2Merge-1;
                end
                %/
                scanStartingSnap=initScanStartingSnap(1:nbrScans2Merge:end);
                if nbrScans2Merge ==1
                  %  fprintf('\nWarning: Merge cancelled --> scans are too large to merge!')
                    flag_scan_merge =0;
                else
                  %  fprintf('--> successful!')
                end
                %/
                
                %/ partitionning
                for idBlk=1:length(scanStartingSnap)-1
                    Blocks(idBlk)=sum(snapshotSzs(scanStartingSnap(idBlk):scanStartingSnap(idBlk+1)-1));
                end
                meanBlkSz=floor(mean(Blocks));
                %handling remaining slice
                lastSlice=sum(snapshotSzs(scanStartingSnap(end):end));
                if  lastSlice >= floor((1-tolBlkSzPc) * min(meanBlkSz,minBlkSz))
                    Blocks=[Blocks(:) ; lastSlice];
                else
                    Blocks(end)=Blocks(end)+lastSlice ;
                end
                %/
            else
               % fprintf("--> adopted\n")
                Blocks = initBlocks;
                flag_scan_merge =0;
                BlkSZ = param.size;
            end
            
        end
        
    end
    %check if further blocking is needed
    meanBlkSz = floor(mean(Blocks));
    stdBlkSz  = std(Blocks)   ;
    if ((meanBlkSz < max(maxBlkSz,BlkSZ)*(1+upperTolBlkSz)) && ...
            (stdBlkSz < 0.5*meanBlkSz) ) || flag_scan_merge
        flag_further_blocking =0;
    end
    
    
    if flag_further_blocking
        fprintf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        fprintf('\n!!Warning: block sizes are large !')
        fprintf('\n!!Warning: re-splitting is needed!')
        fprintf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        %/further blocking
        subParam     = param;
        subParam.pos = Blocks;
        subOutput = util_time_based_block_sp_ar(dataSlice,subParam);
        Blocks = subOutput.block(:);
        fprintf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        fprintf('!\nINFO: re-splitting is DONE       !')
        fprintf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        
    end
    
    output.block =[output.block(:) ;  Blocks(:)];
    output.blockNumber(j) = length(Blocks);
    
    if  verbose
        fprintf("\n %d Blocks, size stats: mean %d, std %d, min %d, max %d",length(Blocks),floor(mean(Blocks)),...
            floor(std(Blocks)),floor(min(Blocks)),floor(max(Blocks)));
    end
    
end
if  verbose
    %     figure(1),plot(abs(output.block-param.size)./param.size,'d')
    output.meanBlkSz = floor(mean(output.block));
    fprintf("\nAverage block size: %d \nDesired block size: %d",output.meanBlkSz,param.size )
    fprintf("\nSmallest block size: %d",min(output.block(:)))
    fprintf("\nLargest block size: %d\n",max(output.block(:)))
end

