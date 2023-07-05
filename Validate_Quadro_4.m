function Validate_Quadro_4(start_day, end_day,excl_days)
% This function validates the data within the files between start_day and end_day, both in terms of missing data, 
% and in the range of the measured variables. 
%
% Version 2 accomodates the three additional houses
%
% Version 3 accomodates the fourth smart plug
% 
% Version 4 accomodates the new inverter variables
%
% If the number of missing values is less than 7, the values are interpolated; if not
% they are left as 0 and the period with no data is marked
%
%   Data is also validated. At present only the ranges of temperature, humidity and
%   solar radiation are verfified. Valid ranges:
%       Smart Pugs: Current [0 inf]
%       WS: AT [-10 50]; RH [0 120]; RAD [0 1500]
%       SPWS Hall: AT [-10 50]; RH [0 120]
%       SPWS Bed 1_2: AT [-10 50]; RH [0 120]
%       SPWS Bed 1_4: AT [-10 50]; RH [0 120]; M [0 100]
%       SPWS L: AT [-10 50]; RH [0 120]; M [0 100]
%       AC: AC_RT [-10 50]; AC_IT [-10 50]; AC_OT [-10 50]
%       
%   Input format
%
%   days should be imput in the format yyyy_MM_dd_HH_mm_ss
%   excl_days is a vector with as many rows as the days excluded
%
%   Output: There is no output from this function. Instead, all data is stored
%   in a file with the same name as the original, but with the last
%   characters 'cor'.
%   
%   Codification of the devices/appliance number:
%       1/1-19: wibees
%       2/20: EM - Original house
%       3/21: Inverter
%       4/22-24: First three Smart Plugs
%       5/25: WS
%       6/26: SPWS Hall
%       6/27: SPWS B12
%       6/28: SPWS B14
%       6/29: SPWS Lounge
%       7/30: AC
%       2/31:EM - New 3P house
%       8/32-33 - New 1P houses
%       4/34 - 4th smart plug
%
%   gap - array with the gaps of the system:record (devices/appliance
%   number/initial sample index;time of the gap beginning/time of the gap
%   end
%   intrec - array with the information of the interpolated data:
%   record(pointer: array with the pointer for the interpolated data for each
%   appliance/devices/appliance/time wehere interpolation starts/time
%   where interpolation ends/duration of the gap (in secs) / number of
%   interpolated vaues
%   intdata - interpolated data: record(k - sample index/t - instant of
%   time/ every relevant interpolated value for the correponding device
%   fault - array with the information of each fault:
%   record(devices/appliance number/relevant variable for the
%   appliance/index of the beginning of the fault/time of the beginning of
%   the fault/ index of the end of the fault/time of the end of
%   the fault
%
%   needed functions: extract_quadro_10


EMF = 0;

if ~isempty(excl_days)
    [m,n]=size(excl_days);
    ed=[];
    for i=1:m-1
        ed=[ed excl_days(i,1:n) '_'];
    end
    ed=[ed excl_days(m,1:n)];
else
    ed='';
end

D=dir;

nfl=length(D);

if ~isempty(excl_days)
    e_d=string(excl_days);
else
    e_d="";
end


nd=0;
for i=1:nfl
    if D(i,1).name~='.'
        if D(i).isdir
            if sum(string(D(i).name(1:10))==e_d)==0
                nd=nd+1;
                dirname(nd,:)=D(i).name;
                dirtime(nd)=datetime(dirname(nd,:),'Format','yyyy_MM_dd_HH_mm_ss');
            end
        end
    end
end

startime=datetime(start_day,'Format','yyyy_MM_dd_HH_mm_ss');
endtime=datetime(end_day,'Format','yyyy_MM_dd_HH_mm_ss');

pos=find(dirtime>=startime);
filestart=pos(1);
if filestart<1
    filestart=1;
end

pos=find(dirtime<=endtime);
filend=max(pos);
clear EMF

for f=filestart:filend
    f1=f-filestart+1;
    disp(' ')
    disp(dirname(f,:))
    disp(' ')
    eval(['cd ',dirname(f,:)])
    ff=[dirname(f,:),'.mat'];
    if exist(ff,'file')
        eval(['load ',ff]);
    else
        cd ..
        extract_quadro_10(dirname(f,:))
        eval(['cd ',dirname(f,:)])
        eval(['load ',ff]);
    end
   
    
    % find missing data
    nEM=0;
    nEM1P=0;
    nINV=0;
    new=0;
    newnew=0;
    nSP=0;
    nWS=0;
    nSPWS_hall1=0;
    nSPWS_bedroom_1_2=0;
    nSPWS_bedroom_1_4=0;
    nSPWS_lounge=0;
    nAC=0;
    nrec=0;
    intrec=[];
    intdata=[];
    % nrec is the number of records in intrec
    ngap=0;
    gap=[];
    % ngap is the number of gaps that could not be interpolated - in the gap
    % record
    nfault=0;
    fault=[];
    % nfault is the number of faults found 

    % Wibees
    [nsamples,ndevices]=size(dt);
    intrec.pointer=zeros(ndevices,1);
    vars={'V','I','F','AP','RP','ApP','PF','AE','IRE','CRE'};

    for i=1:ndevices
        if ndt(i)==0
            [gap,ngap] = no_signal(gap,ngap,1,i,startime,endtime);
        else
            diffdt=seconds(diff(dt(1:ndt(i),i)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
            
            signals=[V(:,i) I(:,i) F(:,i) AP(:,i) RP(:,i) ApP(:,i) PF(:,i) AE(:,i) IRE(:,i) CRE(:,i)];
            clear signalscor;
            [ndtcor(i),dtcor(1:ndtcor(i),i),signalscor(1:ndtcor(i),:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,1,i,dt(:,i),ndt(i),vars,signals,nrec,intrec,intdata,ngap,gap);
            Vcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),1);
            Icor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),2);
            Fcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),3);
            APcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),4);
            RPcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),5);
            ApPcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),6);
            PFcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),7);
            AEcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),8);
            IREcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),9);
            CREcor(1:ndtcor(i),i)=signalscor(1:ndtcor(i),10);
        end
    end
    
    % EM340 - 20
    if exist('dtEM','var') && (ndtEM(1)~=0)
        nEM=1;
        nsamples=length(dtEM);
        vars={'EMVL1_L2','EMVL1_N','EMVL2_L3','EMVL2_N','EMVL3_L1','EMVL3_N','EMVL_N_sys','EMVL_L_sys',...
            'EMI_L1','EMI_L2','EMI_L3','EMAP_sys','EMAP_L1','EMAP_L2','EMAP_L3','EMApP_sys','EMApP_L1','EMApP_L2','EMApP_L3',...
            'EMRP_sys','EMRP_L1','EMRP_L2','EMRP_L3','EMPF_sys','EMPF_L1','EMPF_L2','EMPF_L3','EMF','EMAE_L1','EMAE_L2','EMAE_L3',...
            'EMAET','EMAEP','EMRET','EMREP','EMDP','EMDPP'};

        diffdt=seconds(diff(dtEM(1:ndtEM(1),1)));
        meandiffdt=mean(diffdt);
        k=find(diffdt>=2*meandiffdt);
        mk=length(k);

        if ~exist('EMI_L2','var')
            EMI_L2=EMI_L1;
        end
        signals=[EMVL1_L2(:,1) EMVL1_N(:,1) EMVL2_L3(:,1) EMVL2_N(:,1) EMVL3_L1(:,1) EMVL3_N(:,1) EMVL_N_sys(:,1) EMVL_L_sys(:,1),...
        EMI_L1(:,1) EMI_L2(:,1) EMI_L3(:,1) EMAP_sys(:,1) EMAP_L1(:,1) EMAP_L2(:,1) EMAP_L3(:,1) EMApP_sys(:,1) EMApP_L1(:,1) EMApP_L2(:,1) EMApP_L3(:,1),...
        EMRP_sys(:,1) EMRP_L1(:,1) EMRP_L2(:,1) EMRP_L3(:,1) EMPF_sys(:,1) EMPF_L1(:,1) EMPF_L2(:,1) EMPF_L3(:,1) EMF(:,1) EMAE_L1(:,1) EMAE_L2(:,1) EMAE_L3(:,1),...
        EMAET(:,1) EMAEP(:,1) EMRET(:,1) EMREP(:,1) EMDP(:,1) EMDPP(:,1)];


        clear signalscor;
        [ndtEMcor(1),dtEMcor(1:ndtEMcor(1),1),signalscor(1:ndtEMcor(1),:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,2,ndevices+1,dtEM(:,1),ndtEM(1),vars,signals,nrec,intrec,intdata,ngap,gap);
        EMVL1_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),1);
        EMVL1_Ncor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),2);
        EMVL2_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),3);
        EMVL2_Ncor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),4);
        EMVL3_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),5);
        EMVL3_Ncor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),6);
        EMVL_N_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),7);
        EMVL_L_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),8);

        EMI_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),9);
        EMI_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),10);
        EMI_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),11);
        EMAP_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),12);
        EMAP_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),13);
        EMAP_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),14);
        EMAP_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),15);
        EMApP_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),16);
        EMApP_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),17);
        EMApP_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),18);
        EMApP_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),19);

        EMRP_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),20);
        EMRP_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),21);
        EMRP_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),22);
        EMRP_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),23);
        EMPF_syscor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),24);
        EMPF_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),25);
        EMPF_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),26);
        EMPF_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),27);
        EMFcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),28);
        EMAE_L1cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),29);
        EMAE_L2cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),30);
        EMAE_L3cor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),31);

        EMAETcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),32);
        EMAEPcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),33);
        EMRETcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),34);
        EMREPcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),35);
        EMDPcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),36);
        EMDPPcor(1:ndtEMcor(1),1)=signalscor(1:ndtEMcor(1),37);
    else
        [gap,ngap] = no_signal(gap,ngap,2,20,startime,endtime);
    end
    
    % EM340 - 31
    [~,n]=size(dtEM);
    if n==1
        ndtEM(2)=0;
    end
        
    if ndtEM(2)==0
        [gap,ngap] = no_signal(gap,ngap,2,31,startime,endtime);
    else
        if n==2
            nEM=2;
            nsamples=length(dtEM);
            vars={'EMVL1_L2','EMVL1_N','EMVL2_L3','EMVL2_N','EMVL3_L1','EMVL3_N','EMVL_N_sys','EMVL_L_sys',...
                'EMI_L1','EMI_L2','EMI_L3','EMAP_sys','EMAP_L1','EMAP_L2','EMAP_L3','EMApP_sys','EMApP_L1','EMApP_L2','EMApP_L3',...
                'EMRP_sys','EMRP_L1','EMRP_L2','EMRP_L3','EMPF_sys','EMPF_L1','EMPF_L2','EMPF_L3','EMF','EMAE_L1','EMAE_L2','EMAE_L3',...
                'EMAET','EMAEP','EMRET','EMREP','EMDP','EMDPP'};
    
            diffdt=seconds(diff(dtEM(1:ndtEM(2),2)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
    
            if ~exist('EMI_L2','var')
                EMI_L2=EMI_L1;
            end
            signals=[EMVL1_L2(:,2) EMVL1_N(:,2) EMVL2_L3(:,2) EMVL2_N(:,2) EMVL3_L1(:,2) EMVL3_N(:,2) EMVL_N_sys(:,2) EMVL_L_sys(:,2),...
            EMI_L1(:,2) EMI_L2(:,2) EMI_L3(:,2) EMAP_sys(:,2) EMAP_L1(:,2) EMAP_L2(:,2) EMAP_L3(:,2) EMApP_sys(:,2) EMApP_L1(:,2) EMApP_L2(:,2) EMApP_L3(:,2),...
            EMRP_sys(:,2) EMRP_L1(:,2) EMRP_L2(:,2) EMRP_L3(:,2) EMPF_sys(:,2) EMPF_L1(:,2) EMPF_L2(:,2) EMPF_L3(:,2) EMF(:,2) EMAE_L1(:,2) EMAE_L2(:,2) EMAE_L3(:,2),...
            EMAET(:,2) EMAEP(:,2) EMRET(:,2) EMREP(:,2) EMDP(:,2) EMDPP(:,2)];
    
    
            clear signalscor;
            [ndtEMcor(2),dtEMcor(1:ndtEMcor(2),2),signalscor(1:ndtEMcor(2),:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,2,31,dtEM(:,2),ndtEM(2),vars,signals,nrec,intrec,intdata,ngap,gap);
            EMVL1_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),1);
            EMVL1_Ncor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),2);
            EMVL2_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),3);
            EMVL2_Ncor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),4);
            EMVL3_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),5);
            EMVL3_Ncor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),6);
            EMVL_N_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),7);
            EMVL_L_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),8);
    
            EMI_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),9);
            EMI_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),10);
            EMI_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),11);
            EMAP_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),12);
            EMAP_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),13);
            EMAP_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),14);
            EMAP_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),15);
            EMApP_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),16);
            EMApP_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),17);
            EMApP_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),18);
            EMApP_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),19);
    
            EMRP_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),20);
            EMRP_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),21);
            EMRP_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),22);
            EMRP_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),23);
            EMPF_syscor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),24);
            EMPF_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),25);
            EMPF_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),26);
            EMPF_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),27);
            EMFcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),28);
            EMAE_L1cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),29);
            EMAE_L2cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),30);
            EMAE_L3cor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),31);
    
            EMAETcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),32);
            EMAEPcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),33);
            EMRETcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),34);
            EMREPcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),35);
            EMDPcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),36);
            EMDPPcor(1:ndtEMcor(2),2)=signalscor(1:ndtEMcor(2),37);
        end
    end

    % Inverter - 21
    if exist('dtINV')
        if ndtINV==0
             [gap,ngap] = no_signal(gap,ngap,3,21,startime,endtime);
        else
            nsamples=length(dtINV);
            nINV=1;
            vars={'INVPMPF','INVPMF','INVPMI_L1','INVPMI_L2','INVPMI_L3','INVPMAP_sys','INVPMAP_L1','INVPMAP_L2','INVPMAP_L3',...
                'INVDCI_L1','INVDCI_L2','INVPMRP_L2','INVPMRP_L3','INVPMApP_sys','INVPMApP_L1','INVPMApP_L2','INVPMApP_L3','INVPMV_L1','INVPMV_L2','INVPMV_L3',...
                'INVDCI_L1','INVDCI_L2','INVDCI_L3','INVDCV_L1','INVDCV_L2','INVDCV_L3','INVDCP_L1','INVDCP_L2','INVDCP_L3','INVPF','INVF',...
                'INVI_L1','INVI_L2','INVI_L3','INVAP_sys','INVAP_L1','INVAP_L2','INVAP_L3','INVRP_sys','INVApP_sys','INVV_L1','INVV_L2','INVV_L3',...
                'INVBSC','INVBCC','INVBDC'};
    
            diffdt=seconds(diff(dtINV(1:ndtINV)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);   
    
            signals=[INVPMPF INVPMF INVPMI_L1 INVPMI_L2 INVPMI_L3 INVPMAP_sys INVPMAP_L1 INVPMAP_L2 INVPMAP_L3,...
                INVPMRP_sys INVPMRP_L1 INVPMRP_L2 INVPMRP_L3 INVPMApP_sys INVPMApP_L1 INVPMApP_L2 INVPMApP_L3 INVPMV_L1 INVPMV_L2 INVPMV_L3,...
                INVDCI_L1 INVDCI_L2 INVDCI_L3 INVDCV_L1 INVDCV_L2 INVDCV_L3 INVDCP_L1 INVDCP_L2 INVDCP_L3 INVPF INVF,...
                INVI_L1 INVI_L2 INVI_L3 INVAP_sys INVAP_L1 INVAP_L2 INVAP_L3 INVRP_sys INVApP_sys INVV_L1 INVV_L2 INVV_L3,...
                INVBSC INVBCC INVBDC];
    
            if exist('INVDCP_sys')
                new=1;
                vars=[vars 'INVDCP_sys','INVBGC','INVBT','INVBNC','INVBASC','INVBV',...
                    'INVECB','INVECG','INVECPV','INVECT','INVPCB','INVPCG','INVPCPV',...
                    'INVYD','INVYM','INVYY','INVYT','INVMS','INVIS','INVPL','INVWT','INVCR'];
                signals=[signals INVDCP_sys INVBGC INVBT INVBNC INVBASC INVBV,...
                    INVECB INVECG INVECPV INVECT INVPCB INVPCG INVPCPV,...
                    INVYD INVYM INVYY INVYT INVMS INVIS INVPL INVWT INVCR];
            end
    
            if exist('INVGP')
                newnew=1;
                vars=[vars 'INVGP','INVGE','INVBCDP'];
                signals=[signals INVGP INVGE INVBCDP];
            end

            if exist('INVBatAPSP')
                Cont=1;
                vars=[vars 'INVBatAPSP','INVBatRPSP','INVDeltaCos','INVTotDCcharge','INVTotDCdischarge','INVTotACcharge','INVTotACdischarge','INVTotDCPV','INVTotDCPV1',...
                    'INVTotDCPV2','INVTotDCPV3','INVTotACenergy','INVTotDCpower','INVTotACchargegrid'];
                signals=[signals INVBatAPSP INVBatRPSP INVDeltaCos INVTotDCcharge INVTotDCdischarge INVTotACcharge INVTotACdischarge INVTotDCPV INVTotDCPV1 ...
                    INVTotDCPV2 INVTotDCPV3 INVTotACenergy INVTotDCpower INVTotACchargegrid];
            else
                Cont=0;
            end

    
            clear signalscor;
            [ndtINVcor,dtINVcor(1:ndtINVcor),signalscor(1:ndtINVcor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,3,ndevices+2,dtINV(:),ndtINV,vars,signals,nrec,intrec,intdata,ngap,gap);
            INVPMPFcor(1:ndtINVcor)=signalscor(1:ndtINVcor,1);
            INVPMFcor(1:ndtINVcor)=signalscor(1:ndtINVcor,2);
            INVPMI_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,3);
            INVPMI_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,4);
            INVPMI_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,5);
            INVPMAP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,6);
            INVPMAP_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,7);
            INVPMAP_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,8);
            INVPMAP_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,9);
    
            INVPMRP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,10);
            INVPMRP_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,11);
            INVPMRP_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,12);
            INVPMRP_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,13);
            INVPMApP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,14);
            INVPMApP_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,15);
            INVPMApP_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,16);
            INVPMApP_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,17);
            INVPMV_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,18);
            INVPMV_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,19);
            INVPMV_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,20);
    
            INVDCI_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,21);
            INVDCI_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,22);
            INVDCI_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,23);
            INVDCV_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,24);
            INVDCV_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,25);
            INVDCV_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,26);
            INVDCP_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,27);
            INVDCP_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,28);
            INVDCP_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,29);
            INVPFcor(1:ndtINVcor)=signalscor(1:ndtINVcor,30);
            INVFcor(1:ndtINVcor)=signalscor(1:ndtINVcor,31);
    
            INVI_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,32);
            INVI_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,33);
            INVI_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,34);
            INVAP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,35);
            INVAP_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,36);
            INVAP_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,37);
            INVAP_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,38);
            INVRP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,39);
            INVApP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,40);
            INVV_L1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,41);
            INVV_L2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,42);
            INVV_L3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,43);
    
            INVBSCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,44);
            INVBCCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,45);
            INVBDCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,46);
    
            if exist('INVDCP_sys')
                INVDCP_syscor(1:ndtINVcor)=signalscor(1:ndtINVcor,47);
                INVBGCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,48);
                INVBTcor(1:ndtINVcor)=signalscor(1:ndtINVcor,49);
                INVBNCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,50);
                INVBASCcor(1:ndtINVcor)=signalscor(1:ndtINVcor,51);
                INVBVcor(1:ndtINVcor)=signalscor(1:ndtINVcor,52);
                INVECBcor(1:ndtINVcor)=signalscor(1:ndtINVcor,53);
                INVECGcor(1:ndtINVcor)=signalscor(1:ndtINVcor,54);
                INVECPVcor(1:ndtINVcor)=signalscor(1:ndtINVcor,55);
                INVECTcor(1:ndtINVcor)=signalscor(1:ndtINVcor,56);
                INVPCBcor(1:ndtINVcor)=signalscor(1:ndtINVcor,57);
                INVPCGcor(1:ndtINVcor)=signalscor(1:ndtINVcor,58);
                INVPCPVcor(1:ndtINVcor)=signalscor(1:ndtINVcor,59);
                INVYDcor(1:ndtINVcor)=signalscor(1:ndtINVcor,60);
                INVYMcor(1:ndtINVcor)=signalscor(1:ndtINVcor,61);
                INVYYcor(1:ndtINVcor)=signalscor(1:ndtINVcor,62);
                INVYTcor(1:ndtINVcor)=signalscor(1:ndtINVcor,63);
                INVMScor(1:ndtINVcor)=signalscor(1:ndtINVcor,64);
                INVIScor(1:ndtINVcor)=signalscor(1:ndtINVcor,65);
                INVPLcor(1:ndtINVcor)=signalscor(1:ndtINVcor,66);
                INVWTcor(1:ndtINVcor)=signalscor(1:ndtINVcor,67);
                INVCRcor(1:ndtINVcor)=signalscor(1:ndtINVcor,68);
            end
    
            if exist('INVGP')
                 INVGPcor(1:ndtINVcor)=signalscor(1:ndtINVcor,69);
                 INVGEcor(1:ndtINVcor)=signalscor(1:ndtINVcor,70);
                 INVBCDPcor(1:ndtINVcor)=signalscor(1:ndtINVcor,71);
            end

            if exist('INVBatAPSP')
                INVBatAPSPcor(1:ndtINVcor)=signalscor(1:ndtINVcor,72);
                INVBatRPSPcor(1:ndtINVcor)=signalscor(1:ndtINVcor,73);
                INVDeltaCoscor(1:ndtINVcor)=signalscor(1:ndtINVcor,74);
                INVTotDCchargecor(1:ndtINVcor)=signalscor(1:ndtINVcor,75);
                INVTotDCdischargecor(1:ndtINVcor)=signalscor(1:ndtINVcor,76);
                INVTotACchargecor(1:ndtINVcor)=signalscor(1:ndtINVcor,77);
                INVTotACdischargecor(1:ndtINVcor)=signalscor(1:ndtINVcor,78);
                INVTotDCPVcor(1:ndtINVcor)=signalscor(1:ndtINVcor,79);
                INVTotDCPV1cor(1:ndtINVcor)=signalscor(1:ndtINVcor,80);
                INVTotDCPV2cor(1:ndtINVcor)=signalscor(1:ndtINVcor,81);
                INVTotDCPV3cor(1:ndtINVcor)=signalscor(1:ndtINVcor,82);
                INVTotACenergycor(1:ndtINVcor)=signalscor(1:ndtINVcor,83);
                INVTotDCpowercor(1:ndtINVcor)=signalscor(1:ndtINVcor,84);
                INVTotACchargegridcor(1:ndtINVcor)=signalscor(1:ndtINVcor,85);

            end
        end
    else
        [gap,ngap] = no_signal(gap,ngap,3,21,startime,endtime);
    end
    
    if exist('dtSP')
        nSP=1;
        % Smartt Plugs - 22 to 24 and 34
        [nsamples,nSPlugs]=size(dtSP);

        vars={'SPV','SPI','SPAP','SPAE','SPRssi','SPOn'};

        for i=1:nSPlugs
            if i==4
                ndev = 34;
            else
                ndev = ndevices+2+i;
            end
            if ndtSP(i)==0
                [gap,ngap] = no_signal(gap,ngap,4,ndev,startime,endtime);
            else
                diffdt=seconds(diff(dtSP(1:ndtSP(i),i)));
                meandiffdt=mean(diffdt);
                k=find(diffdt>=2*meandiffdt);
                if isempty(k)
                    mk = 1;
                else
                    mk=length(k);
                end
                % mk is the number of gaps
    
                signals=[SPV(:,i) SPI(:,i) SPAP(:,i) SPAE(:,i) SPRssi(:,i) SPOn(:,i)];
                clear signalscor;
    
                [ndtSPcor(i),dtSPcor(1:ndtSPcor(i),i),signalscor(1:ndtSPcor(i),:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,4,ndev,dtSP(:,i),ndtSP(i),vars,signals,nrec,intrec,intdata,ngap,gap);
                SPVcor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),1);
                SPIcor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),2);
                SPAPcor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),3);
                SPAEcor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),4);
                SPRssicor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),5);
                SPOncor(1:ndtSPcor(i),i)=signalscor(1:ndtSPcor(i),6);
            end
        end
    end
    
    if exist('dtWS')
    % Weather Station - 25
        nWS=1;

        if ndtWS == 0
             [gap,ngap] = no_signal(gap,ngap,5,25,startime,endtime);
        else
            vars={'WS_RH','WS_AT','WS_RAD'};
    
            diffdt=seconds(diff(dtWS(1:ndtWS)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[WS_RH WS_AT WS_RAD];
            clear signalscor;
            [ndtWScor,dtWScor(1:ndtWScor),signalscor(1:ndtWScor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,5,ndevices+6,dtWS,ndtWS,vars,signals,nrec,intrec,intdata,ngap,gap);
            WS_RHcor(1:ndtWScor)=signalscor(1:ndtWScor,1);
            WS_ATcor(1:ndtWScor)=signalscor(1:ndtWScor,2);
            WS_RADcor(1:ndtWScor)=signalscor(1:ndtWScor,3);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,5,25,startime,endtime);
    end
    
        
    if exist('dtH1')
    % SPWS Hall 1 - 26

        nSPWS_hall1=1;
        vars={'H1_AT','H1_RH','H1_L'};
        if ndtH1 == 0
             [gap,ngap] = no_signal(gap,ngap,6,26,startime,endtime);
        else

            diffdt=seconds(diff(dtH1(1:ndtH1)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[H1_AT H1_RH H1_L];
            clear signalscor;
            [ndtH1cor,dtH1cor(1:ndtH1cor),signalscor(1:ndtH1cor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,6,ndevices+7,dtH1,ndtH1,vars,signals,nrec,intrec,intdata,ngap,gap);
            H1_ATcor(1:ndtH1cor)=signalscor(1:ndtH1cor,1);
            H1_RHcor(1:ndtH1cor)=signalscor(1:ndtH1cor,2);
            H1_Lcor(1:ndtH1cor)=signalscor(1:ndtH1cor,3);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,6,26,startime,endtime);
    end
    
    if exist('dtB12')
    % SPWS B12 - 27
        nSPWS_bedroom_1_2=1;
        vars={'B12_AT','B12_RH','B12_L'};
        if ndtB12 == 0
             [gap,ngap] = no_signal(gap,ngap,6,27,startime,endtime);
        else
            diffdt=seconds(diff(dtB12(1:ndtB12)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[B12_AT B12_RH B12_L];
            clear signalscor;
            [ndtB12cor,dtB12cor(1:ndtB12cor),signalscor(1:ndtB12cor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,6,ndevices+8,dtB12,ndtB12,vars,signals,nrec,intrec,intdata,ngap,gap);
            B12_ATcor(1:ndtB12cor)=signalscor(1:ndtB12cor,1);
            B12_RHcor(1:ndtB12cor)=signalscor(1:ndtB12cor,2);
            B12_Lcor(1:ndtB12cor)=signalscor(1:ndtB12cor,3);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,6,27,startime,endtime);
    end
    
    if exist('dtB14')
    % SPWS B14 - 28

        nSPWS_bedroom_1_4=1;
        vars={'B14_AT','B14_RH','B14_WT','B14_M'};

        if ndtB14 == 0
             [gap,ngap] = no_signal(gap,ngap,6,28,startime,endtime);
        else

            diffdt=seconds(diff(dtB14(1:ndtB14)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[B14_AT B14_RH B14_WT B14_M];
            clear signalscor;
            [ndtB14cor,dtB14cor(1:ndtB14cor),signalscor(1:ndtB14cor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,6,ndevices+9,dtB14,ndtB14,vars,signals,nrec,intrec,intdata,ngap,gap);
            B14_ATcor(1:ndtB14cor)=signalscor(1:ndtB14cor,1);
            B14_RHcor(1:ndtB14cor)=signalscor(1:ndtB14cor,2);
            B14_WTcor(1:ndtB14cor)=signalscor(1:ndtB14cor,3);
            B14_Mcor(1:ndtB14cor)=signalscor(1:ndtB14cor,4);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,6,28,startime,endtime);
    end
    
    if exist('dtL')
    % SPWS L - 29
        nSPWS_lounge=1;
        vars={'L_AT','L_RH','L_M'};
        if ndtL == 0
             [gap,ngap] = no_signal(gap,ngap,6,29,startime,endtime);
        else

            diffdt=seconds(diff(dtL(1:ndtL)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[L_AT L_RH L_M];
            clear signalscor;
            [ndtLcor,dtLcor(1:ndtLcor),signalscor(1:ndtLcor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,6,ndevices+10,dtL,ndtL,vars,signals,nrec,intrec,intdata,ngap,gap);
            L_ATcor(1:ndtLcor)=signalscor(1:ndtLcor,1);
            L_RHcor(1:ndtLcor)=signalscor(1:ndtLcor,2);
            L_Mcor(1:ndtLcor)=signalscor(1:ndtLcor,3);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,6,29,startime,endtime);
    end
    
    if exist('dtAC')
    % AC - 30

        nAC=1;
        vars={'AC_PS','AC_SM','AC_EM','AC_TM','AC_OM','AC_RT','AC_IT','AC_OT','AC_FS'};

        if ndtAC == 0
             [gap,ngap] = no_signal(gap,ngap,7,30,startime,endtime);
        else

            diffdt=seconds(diff(dtAC(1:ndtAC)));
            meandiffdt=mean(diffdt);
            k=find(diffdt>=2*meandiffdt);
            mk=length(k);
            % mk is the number of gaps
    
            signals=[AC_PS AC_SM AC_EM AC_TM AC_OM AC_RT AC_IT AC_OT AC_FS];
            clear signalscor;
            [ndtACcor,dtACcor(1:ndtACcor),signalscor(1:ndtACcor,:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,7,ndevices+11,dtAC,ndtAC,vars,signals,nrec,intrec,intdata,ngap,gap);
            AC_PScor(1:ndtACcor)=signalscor(1:ndtACcor,1);
            AC_SMcor(1:ndtACcor)=signalscor(1:ndtACcor,2);
            AC_EMcor(1:ndtACcor)=signalscor(1:ndtACcor,3);
            AC_TMcor(1:ndtACcor)=signalscor(1:ndtACcor,4);
            AC_OMcor(1:ndtACcor)=signalscor(1:ndtACcor,5);
            AC_RTcor(1:ndtACcor)=signalscor(1:ndtACcor,6);
            AC_ITcor(1:ndtACcor)=signalscor(1:ndtACcor,7);
            AC_OTcor(1:ndtACcor)=signalscor(1:ndtACcor,8);
            AC_FScor(1:ndtACcor)=signalscor(1:ndtACcor,9);
        end
    else
        [gap,ngap] = no_signal(gap,ngap,7,30,startime,endtime);
    end

    % EM3401P
    if exist('dtEM1P','var')
        [nsamples,ndevices]=size(dtEM1P);
        vars={'EMV','EMI','EMAP','EMApP','EMRP','EMPF','EMF1P','EMAE1P','EMAEP1P','EMRE1P','EMREP1P','EMDP1P','EMDPP1P'};
        nEM1P=ndevices;
    
        for i=1:ndevices
            if ndtEM1P(i) == 0
                [gap,ngap] = no_signal(gap,ngap,8,31+i,startime,endtime);
            else
                diffdt=seconds(diff(dtEM1P(1:ndtEM1P(i),i)));
                meandiffdt=mean(diffdt);
                k=find(diffdt>=2*meandiffdt);
                mk=length(k);
                % mk is the number of gaps
                
                signals=[EMV(:,i) EMI(:,i) EMAP(:,i) EMApP(:,i) EMRP(:,i) EMPF(:,i) EMF1P(:,i) EMAE1P(:,i) EMAEP1P(:,i) EMRE1P(:,i) EMREP1P(:,i) EMDP1P(:,i) EMDPP1P(:,i)];
                clear signalscor;
                [ndtEM1Pcor(i),dtEM1Pcor(1:ndtEM1Pcor(i),i),signalscor(1:ndtEM1Pcor(i),:),nrec,intrec,intdata,ngap,gap]=update(diffdt,meandiffdt,k,mk,8,31+i,dtEM1P(:,i),ndtEM1P(i),vars,signals,nrec,intrec,intdata,ngap,gap);
                if ndtEM1Pcor(i)==0
                    maxsample=max(ndtcor);
                    pos=find(ndtcor==maxsample);
                    ngap=ngap+1;
                    gap(ngap).devices=8;
                    gap(ngap).num=31+i;
                    gap(ngap).k=1;
                    gap(ngap).tbeg=dtcor(1,1);
                    gap(ngap).tend=dtcor(maxsample,pos(1));
                else
                    EMVcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),1);
                    EMIcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),2);
                    EMAPcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),3);
                    EMApPcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),4);
                    EMRPcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),5);
                    EMPFcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),6);
                    EMF1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),7);
                    EMAE1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),8);
                    EMAEP1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),9);
                    EMRE1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),10);
                    EMREP1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),11);
                    EMDP1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),12);
                    EMDPP1Pcor(1:ndtEM1Pcor(i),i)=signalscor(1:ndtEM1Pcor(i),13);
                end
            end
        end
    else
        [gap,ngap] = no_signal(gap,ngap,8,32,startime,endtime);
        [gap,ngap] = no_signal(gap,ngap,8,33,startime,endtime);
    end
    
    % Validations
    
    %Smart Plugs
    if nSP
        for i=1:nSPlugs
            if i<4
                if ndtSP(i)>0
                    fault=detfault(fault,SPAPcor(:,i),0,+inf, 4,21+i,dtSPcor(:,i),ndtSPcor(i),'SPAP');
                end
            else
                if ndtSP(4)>0
                    fault=detfault(fault,SPAPcor(:,i),0,+inf, 4,34,dtSPcor(:,i),ndtSPcor(i),'SPAP');
                end
            end
        end
    end
    
    %Weather Station
    if nWS && ndtWS>0
        fault=detfault(fault,WS_ATcor,-10,+50, 5,25,dtWScor,ndtWScor,'WS_AT');
        fault=detfault(fault,WS_RHcor,0,+120, 5,25,dtWScor,ndtWScor,'WS_RH');
        fault=detfault(fault,WS_RADcor,0,1500, 5,25,dtWScor,ndtWScor,'WS_RAD');
    end
    
    %SPWS
    if nSPWS_hall1 && ndtH1>0
        fault=detfault(fault,H1_ATcor,-10,+50, 6,26,dtH1cor,ndtH1cor,'H1_AT');
        fault=detfault(fault,H1_RHcor,0,+120, 6,26,dtH1cor,ndtH1cor,'H1_RH');
    end

    if nSPWS_bedroom_1_2 && ndtB12>0
        fault=detfault(fault,B12_ATcor,-10,+50, 6,27,dtB12cor,ndtB12cor,'B12_AT');
        fault=detfault(fault,B12_RHcor,0,+120, 6,27,dtB12cor,ndtB12cor,'B12_RH');
    end            

    if nSPWS_bedroom_1_4 && ndtB14>0
        fault=detfault(fault,B14_ATcor,-10,+50, 6,28,dtB14cor,ndtB14cor,'B14_AT');
        fault=detfault(fault,B14_RHcor,0,+120, 6,28,dtB14cor,ndtB14cor,'B14_RH');
        fault=detfault(fault,B14_Mcor,0,+100, 6,28,dtB14cor,ndtB14cor,'B14_M');
    end     
    
    if nSPWS_lounge && ndtL>0
        fault=detfault(fault,L_ATcor,-10,+50, 6,29,dtLcor,ndtLcor,'L_AT');
        fault=detfault(fault,L_RHcor,0,+120, 6,29,dtLcor,ndtLcor,'L_RH');
        fault=detfault(fault,L_Mcor,0,+100, 6,29,dtLcor,ndtLcor,'L_M');
    end 
    
    if nAC && ndtAC>0
        fault=detfault(fault,AC_RTcor,-10,+50, 7,30,dtACcor,ndtACcor,'AC_RT');
        fault=detfault(fault,AC_ITcor,-10,+50, 7,30,dtACcor,ndtACcor,'AC_IT');
        fault=detfault(fault,AC_OTcor,-10,+50, 7,30,dtACcor,ndtACcor,'AC_OT');
    end
    
    pff=strfind(ff,'.mat');
    lff=length(ff);
    ffcor=[ff(1:lff-4),'_cor.mat'];

    save(ffcor,'ndtcor','dtcor','Vcor','Icor','Fcor','APcor','RPcor','ApPcor','PFcor','AEcor','IREcor','CREcor');
    if nEM>0
        save(ffcor,'EMVL1_L2cor', 'EMVL1_Ncor','EMVL2_L3cor', 'EMVL2_Ncor', 'EMVL3_L1cor', 'EMVL3_Ncor', 'EMVL_N_syscor', 'EMVL_L_syscor', 'EMI_L1cor', 'EMI_L3cor', 'EMAP_syscor', 'EMAP_L1cor', 'EMAP_L2cor', ...
        'EMAP_L3cor', 'EMApP_syscor', 'EMApP_L1cor', 'EMApP_L2cor', 'EMApP_L3cor', 'EMRP_syscor', 'EMRP_L1cor', 'EMRP_L2cor', 'EMRP_L3cor', 'EMPF_syscor', 'EMPF_L1cor', 'EMPF_L2cor', 'EMPF_L3cor', ...
        'EMFcor', 'EMAE_L1cor', 'EMAE_L2cor', 'EMAE_L3cor', 'EMAETcor', 'EMAEPcor', 'EMRETcor', 'EMREPcor', 'EMDPcor', 'EMDPPcor', 'dtEMcor', 'ndtEMcor', '-append');
        if exist('EMI_L2cor','var')
            save(ffcor,'EMI_L2cor','-append')
        end
    end
    if nEM1P>0
        save(ffcor,'EMVcor','EMIcor','EMAPcor','EMApPcor','EMRPcor','EMPFcor','EMF1Pcor','EMAE1Pcor','EMAEP1Pcor','EMRE1Pcor','EMREP1Pcor','EMDP1Pcor','EMDPP1Pcor','dtEM1Pcor','ndtEM1Pcor','-append')
    end
    if nINV>0
        save(ffcor,'INVPMPFcor','INVPMFcor','INVPMI_L1cor','INVPMI_L2cor','INVPMI_L3cor','INVPMAP_syscor','INVPMAP_L1cor','INVPMAP_L2cor','INVPMAP_L3cor','INVPMRP_syscor','INVPMRP_L1cor','INVPMRP_L2cor','INVPMRP_L3cor','INVPMApP_syscor','INVPMApP_L1cor','INVPMApP_L2cor','INVPMApP_L3cor',...
            'INVPMV_L1cor','INVPMV_L2cor','INVPMV_L3cor','INVDCI_L1cor','INVDCI_L2cor','INVDCI_L3cor','INVDCV_L1cor', 'INVDCV_L2cor','INVDCV_L3cor','INVDCP_L1cor','INVDCP_L2cor','INVDCP_L3cor','INVPFcor', 'INVFcor', 'INVI_L1cor', 'INVI_L2cor', 'INVI_L3cor', 'INVAP_syscor', 'INVAP_L1cor', ...
            'INVAP_L2cor', 'INVAP_L3cor', 'INVRP_syscor', 'INVApP_syscor', 'INVV_L1cor', 'INVV_L2cor', 'INVV_L3cor', 'INVBSCcor','INVBCCcor', 'INVBDCcor', 'dtINVcor', 'ndtINVcor', '-append');
        if new
            save(ffcor, 'INVDCP_syscor', 'INVBGCcor', 'INVBTcor', 'INVBNCcor', 'INVBVcor', 'INVECBcor', 'INVECGcor', 'INVECPVcor', 'INVECTcor', 'INVPCBcor', 'INVPCGcor', 'INVPCPVcor', 'INVYDcor', 'INVYMcor', 'INVYYcor', 'INVYTcor', 'INVMScor', 'INVIScor', 'INVPLcor', 'INVWTcor', 'INVCRcor','-append')
            if exist('INVBASCcor','var')
                save(ffcor,'INVBASCcor','-append')
            end
            if newnew
                save(ffcor, 'INVGPcor','INVGEcor', 'INVBCDPcor', '-append')
            end
            if Cont
                save(ffcor, 'INVBatAPSPcor','INVBatRPSPcor','INVDeltaCoscor','INVTotDCchargecor','INVTotDCdischargecor','INVTotACchargecor','INVTotACdischargecor',...
                    'INVTotDCPVcor','INVTotDCPV1cor','INVTotDCPV2cor','INVTotDCPV3cor','INVTotACenergycor','INVTotDCpowercor','INVTotACchargegridcor', '-append')
            end

        end
    end

    if nSP>0
        save(ffcor,'SPVcor', 'SPIcor', 'SPAPcor', 'SPAEcor', 'SPRssicor', 'SPOncor', 'dtSPcor', 'ndtSPcor', '-append');
    end

    if nWS>0
        save(ffcor,'WS_RADcor', 'WS_ATcor', 'WS_RHcor', 'dtWScor', 'ndtWScor', '-append');
    end

    if nSPWS_hall1>0
        save(ffcor,'H1_Lcor', 'H1_ATcor', 'H1_RHcor', 'dtH1cor', 'ndtH1cor', '-append');
    end

    if nSPWS_bedroom_1_2>0
        save(ffcor,'B12_Lcor', 'B12_ATcor', 'B12_RHcor', 'dtB12cor', 'ndtB12cor', '-append');
    end

    if nSPWS_bedroom_1_4>0
        save(ffcor,'B14_Mcor', 'B14_WTcor','B14_ATcor', 'B14_RHcor', 'dtB14cor', 'ndtB14cor', '-append');
    end

    if nSPWS_lounge>0
        save(ffcor,'L_Mcor','L_ATcor', 'L_RHcor', 'dtLcor', 'ndtLcor', '-append');
    end

    if nAC>0
        save(ffcor,'AC_PScor','AC_SMcor', 'AC_EMcor', 'AC_TMcor', 'AC_OMcor', 'AC_RTcor', 'AC_ITcor', 'AC_OTcor', 'AC_FScor', 'dtACcor', 'ndtACcor', '-append');
    end
    
    save(ffcor,'gap','intrec','intdata','fault','-append')
    
    cd ..
end


    
function [ndtcor,dtcor,signalcors,nreccor,intreccor,intdatacor,ngapcor,gapcor]=update(diffdt,meandiffdt,k,mk,devices,num,dt,ndt,signames,signals,nrec,intrec,intdata,ngap,gap)

pos=1;
% pos is the pointer for the original time series
poscor=1;
% poscor is the pointer for the corrected time series
[~,nvars]=size(signals);

nrecj=0;
totsamp=0;
ndtcor=ndt;
if (mk>0) && (nrec>1)
    npoint=length(intrec(nrec).pointer);
    nad=num-npoint;
    if num>npoint
        for i=1:nrec
            intrec(i).pointer(npoint+1:num)=zeros(nad,1);
        end
    end
end
        
for j=1:mk  
    dur=diffdt(k(j));
    newsamples=floor(dur/meandiffdt-1);
    % newsamples is the number of samples to interpolate
    if newsamples<7
        %interpolate
        totsamp=totsamp+newsamples;
        %previous time slice
        poscorn=poscor+k(j)-pos;
        dtcor(poscor:poscorn)=dt(pos:k(j));
        signalcors(poscor:poscorn,:)=signals(pos:k(j),:);

        pos=k(j)+1;
        poscor=poscorn;

        if k(j)<=10
            ini=1;
        else
            ini=k(j)-10;
        end
        if (k(j)+newsamples)>(ndt-10)
            fim=ndt;
        else
            fim=k(j)+newsamples+10;
        end

        nrec=nrec+1;

        intrec(nrec).devices=devices;
        intrec(nrec).num=num;
        intrec(nrec).startgap=dt(k(j));
        intrec(nrec).endgap=dt(k(j)+1);
        intrec(nrec).durgap=diffdt(k(j));
        intrec(nrec).newsamples=newsamples;
        if nrec>1
            intrec(nrec).pointer=intrec(nrec-1).pointer;
            nd=length(intrec(nrec).pointer);
            if nd>=num
                intrec(nrec).pointer(num)=intrec(nrec-1).pointer(num)+intrec(nrec-1).newsamples; 
            else
                intrec(nrec).pointer(num)=intrec(nrec-1).newsamples;
            end
        end
        if nrecj==0
            intrec(nrec).pointer(num)=1;
            nrecj=1;
        end
        pnt=intrec(nrec).pointer(num)-1;

        for ind=1:newsamples                    
            intdata(pnt+ind,num).k=poscor+ind;
            intdata(pnt+ind,num).t=dt(k(j))+seconds(meandiffdt)*ind;
            dtcor(poscor+ind)=intdata(pnt+ind,num).t;
        end
        ndtcor=ndtcor+newsamples;

        %   interpolate
        for i=1:nvars
            y=[signals(ini:k(j),i);NaN*ones(newsamples,1);signals(k(j)+1:fim,i)];
            ycleaned=fillmissing(y,'movmedian',[2 2]); 
            for ind=1:newsamples
                eval(['intdata(pnt+ind,num).', signames{i}, '=ycleaned(k(j)-ini+ind+1);'])
                eval(['signalcors(poscor+ind,i)=intdata(pnt+ind,num).', signames{i},';'])
            end
        end
        
        poscor=poscor+newsamples+1;
    else
        % the gap is too large. 
        ngap=ngap+1;
        gap(ngap).devices=devices;
        gap(ngap).num=num;
        gap(ngap).k=k(j)+totsamp;
        gap(ngap).tbeg=dt(k(j));
        gap(ngap).tend=dt(k(j)+1);
    end
end

if mk>0
    % update the rest of the time series
    dtcor(poscor:ndtcor)=dt(pos:ndt);
    signalcors(poscor:ndtcor,:)=signals(pos:ndt,:);
else
    ndtcor=ndt;
    dtcor=dt(1:ndtcor);
    signalcors=signals(1:ndtcor,:);
end

dtcor=dtcor(:);
nreccor=nrec;
intreccor=intrec;
intdatacor=intdata;
ngapcor=ngap;
gapcor=gap;

function fault=detfault(fault,signal,lb,up, devices,num,dt,ndt,var)

nfault=length(fault);

if lb==-inf
    kfault=find(signal>ub);
elseif up==+inf
    kfault=find(signal<lb);
else
    kfault=union(find(signal<lb),find(signal>up));
end
if ~isempty(kfault)
    %there is a fault
    nfault=nfault+1;
    fault(nfault).devices=devices;
    fault(nfault).num=num;
    fault(nfault).var={var};
    fault(nfault).kbeg=kfault(1);
    fault(nfault).tbeg=dt(kfault(1));
    nf=length(kfault);
    fim=kfault(nf)+1;
    if fim+ndt*0.005>ndt
        fim=ndt;
    end
    fault(nfault).kend=fim;
    fault(nfault).tend=dt(fim);
end

function [gap,ngap] = no_signal(gap,ngap,devices, num,startime,endtime)

ngap = ngap + 1;
gap(ngap).devices = devices;
gap(ngap).num = num;
gap(ngap).k = 1;
gap(ngap).tbeg = startime;
gap(ngap).tend = endtime;
    