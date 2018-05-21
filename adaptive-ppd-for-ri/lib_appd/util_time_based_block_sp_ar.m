function out =util_time_based_block_sp_ar(time_,param)
param.pos=[0 param.pos];
out.partition = [] ;
flag_repeat_blocking =0;
for j =1:length(param.pos)-1
    flag_forced_split=1;
    blockSZ = max(1, floor(param.pos(j+1)/param.size));
    blockMod = mod(param.pos(j+1),param.size);
    if (blockMod > (param.size *0.5)) && (param.size<param.pos(j+1))
        blockSZ = blockSZ+1;
    end    
    last = sum(param.pos(1:j+1));
    first =1;
    if j>1
        first = sum(param.pos(1:j))+1;
    end
    timec = time_(first:last);    
    M = length(timec);
    %%
    t =1;
    cmt = 1;
    snap_sz =[];
    snaps = [];
    start = 2;
    for  i = start: length(timec)
        if timec(i) ~=timec(i-1) %new snapshot
            snap_sz(t)= i-cmt;
            cmt =  i;
            snaps(i-1) = t;
            t=t+1;
        else
            snaps(i-1) =t;
        end
    end
    snaps=[snaps snaps(end)];
    snaps(1:start-1) = snaps(start);
    snap_sz(t) = M - cmt+1;
    
    if length(snap_sz)~=snaps(end)
        error('! CHECK SNAPSHOT DELIMITING')
    end
    %%
    diffSnapTime = [];
    for i=1:length(snap_sz)-1
        
        snapEndTime=timec(sum(snap_sz(1:i)));
        snapStartTime=timec(sum(snap_sz(1:i))+1);
        diffSnapTime(i)  = snapStartTime -snapEndTime;
    end
    [~,startTrack]=find(abs(diffSnapTime)>3*median(abs(diffSnapTime)));
    if isempty(startTrack)
        dummy=  floor(length(snap_sz)/blockSZ);
        startTrack=dummy:dummy:(length(snap_sz));
        startTrack(end)=startTrack(end)+mod(length(snap_sz),blockSZ);
        flag_forced_split =0;
        max_nblocks = blockSZ;
    else
        max_nblocks = length(startTrack)+2;
    end
    
    %%
    slice=0;
    partition =[];
    blockSZ = min(max_nblocks,blockSZ);    
    
    if max_nblocks==1 || blockSZ==1
        partition =sum(snap_sz);       
    else
        t =1;
        i =1;
        snapStart=1;        
        snapLast=startTrack(t);
        if (blockSZ <=max_nblocks) && flag_forced_split==0
            while (t<max_nblocks+1) &&(snapLast-1<length(snap_sz)) && i<1+max_nblocks
                old_slice=slice;
                slice=sum(snap_sz(snapStart:snapLast));
                
                if slice>=1.2*param.size
                    slice = old_slice;
                    i=i+1;
                    partition=[partition slice];
                    snapStart=startTrack(t-1)+1;
                    snapLast=startTrack(t);                  
                    
                elseif slice>param.size && slice<1.2*param.size
                    i=i+1;
                    partition=[partition slice];                    
                    snapStart=snapLast+1;
                    if t> length(startTrack)-1
                        snapLast= length(snap_sz);
                    elseif flag_forced_split
                        snapLast=startTrack(t+1);
                        t=t+1;
                    end
                    
                elseif snapLast ==length(snap_sz)
                    if slice>param.size *0.5
                        partition=[partition slice];
                    else
                        partition(end)=partition(end)+slice;
                    end
                   break;
                    
                elseif  (t<length(startTrack))
                    snapLast= startTrack(t+1);
                    t=t+1;
                    
                elseif (t==length(startTrack))
                    snapLast=length(snap_sz);
                end                
            end           
            
        else            
            partition(1)=sum(snap_sz(1:startTrack(1)));            
            for i=2:length(startTrack)
                partition(i)=sum(snap_sz(startTrack(i-1)+1:startTrack(i)));
            end
            partition(length(startTrack)+1)=sum(snap_sz(startTrack(end)+1:end));
                      
        end
    end
    for i=1:length(partition)
        if partition(i) >2*param.size
            flag_repeat_blocking=1;
        end
    end
    if sum(partition)~=param.pos(j+1)
        if param.pos(j+1)-sum(partition)>param.size*0.5
            partition=[partition  param.pos(j+1)-sum(partition)];
        else
            partition(end)=partition(end)+param.pos(j+1)-sum(partition)  ;
        end        
         error('! CHECK BLOCK SPLITTING')
    end
    if flag_repeat_blocking       
        param.pos=partition;
        outn =util_time_based_block_sp_ar(time_,param);
        out.partition =[out.partition outn.partition];
    else
        out.blockNumber(j) = length(partition);
        out.partition =[out.partition partition];
    end
    
end


end
