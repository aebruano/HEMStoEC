function extract_quadro_10(dir_name, map)
% This function, which should be executed at the directory above the daily measurements, extracts the values measured at each circuit, and groups them by variable. 
% At present, there are 19 separate circuits (16 monophasic + 1 triphasic),
% 1 EM340, 3 Smart Plugs, and one inverter (Plenticore)
% All extracted data are stored in a matlab data file whose name coincides with the directory name.
%
%Version 3 multiplies the current, power and energy values of the wibbes by
%the corresponding factor
%
%Version 4 includes NewVarsVersion
%
%Version 5 includes the Weather Station ans the Self-Powered Wireless
%Sensors
%
%Version 6 optimizes the code and includes the SPWS
%
%Version 7 includes the air conditioner
%
%Version 8 includes the other houses
%
%Version 9 incorporates the 4th SmartPlug
%
%Version 10 incorporates the variables for External Battery Management

% Functions and data files used:
%   Factor.mat

%   Variables names:

%   For Weebeees:

%   V - Voltage
%   I - Current
%   F - Frequency
%   AP - ActivePower
%   RP - Reactive Power
%   ApP - Apparent Power
%   PF - Power Factor
%   AE - Active Energy
%   IRE - Inductive Reactive Energy
%   CRE - Capacitive Reactive Energy
%
%   Only 3P:
%
%   average_F - Avergafe Frequency
%   T_AP - Total Active Power
%   T_RP - Total Reactive Power
%   T_ApP - Total Apparent Power
%   T_PF - Total Power Factor
%   T_AE - Total Active Energy
%   T_IRE - Total Inductive Reactive Energy
%   T_CRE - Total Capacitive Reactive Energy

%   EM340 and EM340_C2: 3P Total House Consumption (1 - minha; 2 - Carlos)
%
%   Voltages:  	
%       EMVL1_L2
%       EMVL1_N
%       EMVL2_L3
%       EMVL2_N
%       EMVL3_L1
%       EMVL3_N
%       EMVL_N_sys
%       EMVL_L_sys
%   Current:		
%       EMI_L1
%       EMI_L2
%       EMI_L3
%   Active Power:	
%       EMAP_sys
%       EMAP_L1
%       EMAP_L2
%       EMAP_L3
%   Apparent Power:	
%       EMApP_sys
%       EMApP_L1
%       EMApP_L2
%       EMApP_L3
%   Reactive Power:	
%       EMRP_sys
%       EMRP_L1
%       EMRP_L2
%       EMRP_L3
%   Power Factor:	
%       EMPF_sys
%       EMPF_L1
%       EMPF_L2
%       EMPF_L3     
%   Harmonic Distortion:	
%       EMTHDV_sys
%       EMTHDV_L1
%       EMTHDV_L2
%       EMTHDV_L3
%       EMTHDI_L1
%       EMTHDI_L2
%       EMTHDI_L3
%   Frequency:
%       EMF			
%   Active Energy:		
%       EMAE_L1
%       EMAE_L2
%       EMAE_L3
%   Total Energy and Demand Power:	
%       EMAET - Active Energy Total
%       EMAEP - Active Energy Partial
%       EMRET - Reactive Energy Total
%       EMREP - Reactive Energy Partial
%       EMDP - Demand Power
%       EMDPP - Demand Power Peak
%
%   EM112 C1 and C3: 1P Total House Consumption
%
%   Voltages:  	
%       EMV
%   Current:		
%       EMI
%   Active Power:	
%       EMAP
%   Apparent Power:	
%       EMApP
%   Reactive Power:	
%       EMRP
%   Power Factor:	
%       EMPF 
%   Frequency:
%       EMF1P			
%   Total Energy and Demand Power:	
%       EMAE1P - Active Energy
%       EMAEP1P - Active Energy Partial
%       EMRE1P - Active Energy
%       EMREP1P - Reactive Energy Partial
%       EMDP1P - Power Demand 
%       EMDPP1P - PowerDemand Peak
%
%   Inverter

%   Powermeter - Supplied by the Grid
%       INVPMPF - Power Factor
%       INVPMF - Frequency
%       Current
%           INVPMI_L1
%           INVPMI_L2
%           INVPMI_L3
%       Active Power
%           INVPMAP_sys
%           INVPMAP_L1
%           INVPMAP_L2
%           INVPMAP_L3
%       Reactive Power
%           INVPMRP_sys
%           INVPMRP_L1
%           INVPMRP_L2
%           INVPMRP_L3
%       Apparent Power
%           INVPMApP_sys
%           INVPMApP_L1
%           INVPMApP_L2
%           INVPMApP_L3
%       Voltage
%           INVPMV_L1
%           INVPMV_L2
%           INVPMV_L3
%
%   DC - DC Part of the inverter
%       Current
%           INVDCI_L1
%           INVDCI_L2
%           INVDCI_L3
%       Voltage
%           INVDCV_L1
%           INVDCV_L2
%           INVDCV_L3
%       Power
%           INVDCP_L1
%           INVDCP_L2
%           INVDCP_L3
%          
%   Others - AC part of the inverter (supplied by PV and Battery
%       INVPF - Power Factor
%       INVF - Frequency
%       Current
%           INVI_L1
%           INVI_L2
%           INVI_L3
%       Active Power
%           INVAP_sys
%           INVAP_L1
%           INVAP_L2
%           INVAP_L3
%       Reactive Power
%           INVRP_sys
%       Apparent Power
%           INVApP_sys
%       Voltage
%           INVV_L1
%           INVV_L2
%           INVV_L3
%       Battery
%           INVBSC - Battery State Charge
%           INVBCC - Battery Charge Current
%           INVBDC - Battery Discharge Current
%
%  New version
%       DC
%           INVDCP_sys - Total Power
%       Battery
%           INVBGC - Gross Capacity
%           INVBT - Temperature
%           INVBNC - Number of Cycles
%           INVBASC - Actual State of Charge
%           INVBV - Voltage
%       Home
%           Energy Consumption
%               INVECB - From Battery
%               INVECG - From Grid
%               INVECPV - From PV
%               INVECT - Total Energy Consumption
%           Power Consumption
%               INVPCB - From Battery
%               INVPCG - From Grid
%               INVPCPV - From PV
%           Yield
%               INVYD - Daily
%               INVYM - Monthly
%               INVYY - Yearly
%               INVYT -- Total
%           Others
%               INVMS - Manager State
%               INVIS - Inverter State
%               INVPL - Power limit
%               INVWT - Work Time
%               INVCR - Total Home Consumption Rate
%               
%   NewVarsVersion
%       INVGP - Invertor Generation Power
%       INVGE - Invertor Generation Energy
%       INVBCDP - Battery Charge/Discharge Power
%
%   Battery Management
%       INVBatAPSP - Active Power Setpoint
%       INVBatRPSP - Reactive Power Setpoint
%       INVDeltaCos - Delta cos Setpoint
%       INVTotDCcharge - Total DC Charge Energy (DC-side-to-battery)
%       INVTotDCdischarge - Total DC Discharge Energy (DC-side-from-battery)
%       INVTotACcharge - Total AC Charge Energy (AC-side-to-battery)
%       INVTotACdischarge - Total AC Discharge Energy (AC-side-from-battery)
%       INVTotDCPV - Total DC energy (Sum of all inputs)
%       INVTotDCPV1 - Total DC energy from PV1
%       INVTotDCPV2 - Total DC energy from PV2
%       INVTotDCPV3 - Total DC energy from PV3
%       INVTotACenergy - Total AC energy AC_side to grid
%       INVTotDCpower - Total DC power
%       INVTotACchargegrid - Total AC charge energy grid
%
%   Smart Plugs
%       SPV - Voltage
%       SPI - Current
%       SPAP - Active Power
%       SPAE - Active Energy
%       SPRssi - Signal Power
%       SPOn - On/Off
%
%   Weather Station:
%       WS_RH - Relative Humidity
%       WS_AT - Air Temperature
%       WS_RAD - Solar Radiation
%
%   Self-Powered Wirelees Sensors
%       Hall 1:
%           H1_AT - Air Temperature
%           H1_RH - Relative Humidity
%           H1_L - Light
%       Bedroom 1_2:
%           B12_AT - Air Temperature
%           B12_RH - Relative Humidity
%           B12_L - Light
%       Bedroom 1_4:
%           B14_AT - Air Temperature
%           B14_RH - Relative Humidity
%           B14_WT - Wall Temperature
%           B14_M - Movement
%       Lounge:
%           L_AT - Air Temperature
%           L_RH - Relative Humidity
%           L_M - Movement
%
%   Air conditioner:
%       AC_PS - Power State (On/off)
%       AC_SM - Swing Mode (On/off)
%       AC_EM - Eco Mode (On/off)
%       AC_TM - Turbo Mode (On/off)
%       AC_OM - Operational Mode (Numerical)
%       AC_RT - Reference Temperature (ºC)
%       AC_IT - Indoor Temperature (ºC)
%       AC_OT - OutDoor Temperature (ºC)
%       AC_FS - Fan Speed (%)



load('Factor.mat','Factor')

eval(['cd ',dir_name]);

if exist([dir_name '.mat'],'file')
    delete([dir_name,'.mat'])
end

texto=importdata('VARS');
ntxt=length(texto);

for i=2:ntxt
    Words(i-1,:)=split(texto(i),',');
end

ntxt=ntxt-1;

 for i=1:ntxt
    sensorlist(i,:)=Words{i,1};
end

D=dir;

nfl=length(D);

n3P=0;
n1P=0;
nEM=0;
nEM1P=0;
nINV=0;
nSP=0;
nWS=0;
nSPWS_hall1=0;
nSPWS_bedroom_1_2=0;
nSPWS_bedroom_1_4=0;
nSPWS_lounge=0;
nAC=0;

for i=1:nfl
    Files{i}=D(i).name;
     if ~isempty(strfind(Files{i},'plenticore'))
        nINV=nINV+1;
        IndINV(nINV)=i;
     elseif ~isempty(strfind(Files{i},'smartplug'))
        number=str2num(Files{i}(17));
        if number>nSP
            nSP=number;
        end
        TabSP(number)=1;
        IndSP(number)=i;
     elseif ~isempty(strfind(Files{i},'em340rtu'))
         if nEM==0
             nEM=1;
         end
        IndEM(1)=i;
     elseif ~isempty(strfind(Files{i},'em340C2_rtu'))
        nEM=2;
        IndEM(2)=i;
     elseif ~isempty(strfind(Files{i},'em112C1_rtu'))
         if nEM1P==0
             nEM1P=1;
         end
        IndEM1P(1)=i;
     elseif ~isempty(strfind(Files{i},'em112C3_rtu'))
        nEM1P=2;
        IndEM1P(nEM1P)=i;
    elseif ~isempty(strfind(Files{i},'3p'))|~isempty(strfind(Files{i},'3P'))
        n3P=n3P+1;
        Ind3P(n3P)=i;
     elseif ~isempty(strfind(Files{i},'1p'))|~isempty(strfind(Files{i},'1P'))
        n1P=n1P+1;
        Ind1P(n1P)=i;
     elseif ~isempty(strfind(Files{i},'ws038798'))
         nWS=nWS+1;
         IndWS(nWS)=i;
     elseif ~isempty(strfind(Files{i},'em_thlbs_c'))
         nSPWS_hall1=nSPWS_hall1+1;
         IndSPWS_hall1(nSPWS_hall1)=i;
     elseif ~isempty(strfind(Files{i},'em_thlbs_qp'))
         nSPWS_bedroom_1_2=nSPWS_bedroom_1_2+1;
         IndSPWS_bedroom_1_2(nSPWS_bedroom_1_2)=i;
     elseif ~isempty(strfind(Files{i},'em_thmwtbs_qg'))
         nSPWS_bedroom_1_4=nSPWS_bedroom_1_4+1;
         IndSPWS_bedroom_1_4(nSPWS_bedroom_1_4)=i;
     elseif ~isempty(strfind(Files{i},'em_thmbs_s'))
         nSPWS_lounge=nSPWS_lounge+1;
         IndSPWS_lounge(nSPWS_lounge)=i;
     elseif ~isempty(strfind(Files{i},'mi_air_cond'))
         nAC=nAC+1;
         IndAC(nAC)=i;
    end
end

disp('1 Phase');

nwb=n1P+3;
ndt=zeros(nwb,1);
maxsamplesdt=24*60*60;
dt=datetime(zeros(maxsamplesdt,nwb),'ConvertFrom','datenum');
V=zeros(maxsamplesdt,nwb);
AP=zeros(maxsamplesdt,nwb);
RP=zeros(maxsamplesdt,nwb);
ApP=zeros(maxsamplesdt,nwb);
PF=zeros(maxsamplesdt,nwb);
AE=zeros(maxsamplesdt,nwb);
IRE=zeros(maxsamplesdt,nwb);
CRE=zeros(maxsamplesdt,nwb);
F=zeros(maxsamplesdt,nwb);
I=zeros(maxsamplesdt,nwb);



for ind=1:n1P
    nome=Files{Ind1P(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    if length(nome)==16
        i=str2double(nome(15:16));
    else
        i=str2double(nome(15));
    end
    
    if nargin==2
        i=map(i); 
    end
    
    if i==18
        i=19;
    end
    
    eval(['instant=',nomevar,'.instant;']);
    ndt(i)=length(instant);
    dt(1:ndt(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if ~isempty(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if ~isempty(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'V'))
                        V(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'AP'))
                        AP(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        AP(1:ndt(i),i)=AP(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'RP'))
                        RP(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        RP(1:ndt(i),i)=RP(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'ApP'))
                        ApP(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        ApP(1:ndt(i),i)=ApP(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'PF'))
                        PF(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);    
                    elseif ~isempty(strfind(Words{k,4},'AE'))
                        AE(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        AE(1:ndt(i),i)=AE(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'IRE'))
                        IRE(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        IRE(1:ndt(i),i)=IRE(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'CRE'))
                        CRE(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        CRE(1:ndt(i),i)=CRE(1:ndt(i),i)*Factor(i);
                    elseif ~isempty(strfind(Words{k,4},'F'))
                        F(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'I'))
                        I(1:ndt(i),i)=eval([nomevar,'.',sensors{j}]);
                        I(1:ndt(i),i)=I(1:ndt(i),i)*Factor(i);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

% 3 phase

pos=n1P;
pos3=1;
ndt3=zeros(n3P,1);
dt3=datetime(zeros(maxsamplesdt,n3P),'ConvertFrom','datenum');
average_F=zeros(maxsamplesdt,n3P);
T_AP=zeros(maxsamplesdt,n3P);
T_RP=zeros(maxsamplesdt,n3P);
T_ApP=zeros(maxsamplesdt,n3P);
T_PF=zeros(maxsamplesdt,n3P);
T_AE=zeros(maxsamplesdt,n3P);
T_IRE=zeros(maxsamplesdt,n3P);
T_CRE=zeros(maxsamplesdt,n3P);
disp('3 Phases');

for ind=1:n3P
    nome=Files{Ind3P(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    ndt(pos)=length(instant);
    ndt(pos+1)=ndt(pos);
    ndt(pos+2)=ndt(pos);
    ndt3(pos3)=ndt(pos);
    dt(1:ndt(pos),pos)=datetime(instant);
    dt(1:ndt(pos+1),pos+1)=dt(1:ndt(pos),pos);
    dt(1:ndt(pos+2),pos+2)=dt(1:ndt(pos),pos);
    dt3(1:ndt3(pos3),pos3)=dt(1:ndt(pos),pos);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    
    for j=1:nsen
        if ~isempty(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if ~isempty(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'V'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            V(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            V(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            V(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                        end
                    elseif ~isempty(strfind(Words{k,4},'average_F'))
                        average_F(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'total_AP'))   
                        T_AP(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'T_RP'))   
                        T_RP(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'totalApP'))   
                        T_ApP(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);  
                    elseif ~isempty(strfind(Words{k,4},'T_PF'))   
                        T_PF(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);  
                    elseif ~isempty(strfind(Words{k,4},'T_AE'))   
                        T_AE(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);     
                    elseif ~isempty(strfind(Words{k,4},'T_IRE'))   
                        T_IRE(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);     
                    elseif ~isempty(strfind(Words{k,4},'T_CRE'))   
                        T_CRE(1:ndt3(pos3),pos3)=eval([nomevar,'.',sensors{j}]);  
                    elseif ~isempty(strfind(Words{k,4},'AP'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            AP(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            AP(1:ndt(pos),pos)=AP(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            AP(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            AP(1:ndt(pos+1),pos+1)=AP(1:ndt(pos+1),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            AP(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            AP(1:ndt(pos+2),pos+2)=AP(1:ndt(pos+2),pos+2)*Factor(pos+2);
                        end
                    elseif ~isempty(strfind(Words{k,4},'RP'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            RP(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            RP(1:ndt(pos),pos)=RP(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            RP(1:ndt(pos+1),pos+1)=eval([nomevar,'.',sensors{j}]);
                            RP(1:ndt(pos+1),pos+1)=RP(1:ndt(pos+1),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            RP(1:ndt(pos+2),pos+2)=eval([nomevar,'.',sensors{j}]);
                            RP(1:ndt(pos+2),pos+2)=RP(1:ndt(pos+2),pos+2)*Factor(pos+2);
                        end
                    elseif ~isempty(strfind(Words{k,4},'ApP'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            ApP(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            ApP(1:ndt(pos),pos)=ApP(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            ApP(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            ApP(1:ndt(pos),pos+1)=ApP(1:ndt(pos),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            ApP(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            ApP(1:ndt(pos),pos+2)=ApP(1:ndt(pos),pos+2)*Factor(pos+2);
                        end

                    elseif ~isempty(strfind(Words{k,4},'PF'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            PF(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            PF(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            PF(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                        end
   
                    elseif ~isempty(strfind(Words{k,4},'AE'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            AE(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            AE(1:ndt(pos),pos)=AE(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            AE(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            AE(1:ndt(pos),pos+1)=AE(1:ndt(pos),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            AE(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            AE(1:ndt(pos),pos+2)=AE(1:ndt(pos),pos+2)*Factor(pos+2);
                        end

                    elseif ~isempty(strfind(Words{k,4},'IRE'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            IRE(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            IRE(1:ndt(pos),pos)=IRE(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            IRE(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            IRE(1:ndt(pos),pos+1)=IRE(1:ndt(pos),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            IRE(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            IRE(1:ndt(pos),pos+2)=IRE(1:ndt(pos),pos+2)*Factor(pos+2);
                        end

                    elseif ~isempty(strfind(Words{k,4},'CRE'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            CRE(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            CRE(1:ndt(pos),pos)=CRE(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            CRE(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            CRE(1:ndt(pos),pos+1)=CRE(1:ndt(pos),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            CRE(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            CRE(1:ndt(pos),pos+2)=CRE(1:ndt(pos),pos+2)*Factor(pos+2);
                        end

                    elseif ~isempty(strfind(Words{k,4},'F'))
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            F(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            F(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            F(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                        end

                    elseif ~isempty(strfind(Words{k,4},'I'))
                        % Current
                        if ~isempty(strfind(Words{k,4},'L1'))                           
                            I(1:ndt(pos),pos)=eval([nomevar,'.',sensors{j}]);
                            I(1:ndt(pos),pos)=I(1:ndt(pos),pos)*Factor(pos);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            I(1:ndt(pos),pos+1)=eval([nomevar,'.',sensors{j}]);
                            I(1:ndt(pos),pos+1)=I(1:ndt(pos),pos+1)*Factor(pos+1);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            I(1:ndt(pos),pos+2)=eval([nomevar,'.',sensors{j}]);
                            I(1:ndt(pos),pos+2)=I(1:ndt(pos),pos+2)*Factor(pos+2);
                        end
                    end
                end
                k=k+1;
            end   
        end
    end
    T_AP(1:ndt3(pos3),pos3)=sum(AP(1:ndt3(pos3),pos:pos+2)');
    T_RP(1:ndt3(pos3),pos3)=sum(RP(1:ndt3(pos3),pos:pos+2)');
    T_ApP(1:ndt3(pos3),pos3)=sum(ApP(1:ndt3(pos3),pos:pos+2)');
    T_AE(1:ndt3(pos3),pos3)=sum(AE(1:ndt3(pos3),pos:pos+2)');
    T_IRE(1:ndt3(pos3),pos3)=sum(IRE(1:ndt3(pos3),pos:pos+2)');
    T_CRE(1:ndt3(pos3),pos3)=sum(CRE(1:ndt3(pos3),pos:pos+2)');

    pos=pos+3;
    pos3=pos3+1;
    eval('cd ..')
end

disp('EM 340 1P')
ndtEM1P=zeros(nEM1P,1);
dtEM1P=datetime(zeros(maxsamplesdt,nEM1P),'ConvertFrom','datenum');  	
EMV=zeros(maxsamplesdt,nEM1P);		
EMI=zeros(maxsamplesdt,nEM1P);
EMAP=zeros(maxsamplesdt,nEM1P);
EMApP=zeros(maxsamplesdt,nEM1P);
EMRP=zeros(maxsamplesdt,nEM1P);
EMPF=zeros(maxsamplesdt,nEM1P);
EMF1P=zeros(maxsamplesdt,nEM1P);		
EMAE1P=zeros(maxsamplesdt,nEM1P);
EMRE1P=zeros(maxsamplesdt,nEM1P);
EMAEP1P=zeros(maxsamplesdt,nEM1P);
EMREP1P=zeros(maxsamplesdt,nEM1P);
EMDP1P=zeros(maxsamplesdt,nEM1P); 
EMDPP1P=zeros(maxsamplesdt,nEM1P);

for ind=1:nEM1P
    if IndEM1P(ind)~=0
        nome=Files{IndEM1P(ind)};
        eval(['cd ',nome]);
        d=dir;
        dev=d(3).name(1:2);
        nomef=[dev, '_',nome];
        eval(['load ',nomef]);
        nomevar=['data_',nomef];
        if nome(6)=='C' % asneira do Sergio
            nome(6)='c';
            nomef=[dev, '_',nome];
            nomevar=['data_',nomef];
        end
        
        eval(['instant=',nomevar,'.instant;']);
        ndtEM1P(ind)=length(instant);
    
        dtEM1P(1:ndtEM1P(ind),ind)=datetime(instant);
        
        eval(['sensors=fieldnames(',nomevar,');']);
        nsen=length(sensors);
     
        for j=1:nsen
            if ~isempty(strmatch(sensors(j),'instant'))==0
                k=1;
                while k<=ntxt
                    if ~isempty(strmatch(sensorlist(k,:),sensors{j}))
                        if ~isempty(strfind(Words{k,4},'_V'))              
                            EMV(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);     
                        elseif ~isempty(strfind(Words{k,4},'_I'))
                            % Current
                            EMI(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);  
                        elseif ~isempty(strfind(Words{k,4},'_PF'))
                            % Power Factor
                            EMPF(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_ApP'))
                            % Apparent Power
                            EMApP(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_RP'))
                            % Reactive Power
                            EMRP(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_PD'))
                           % Power Demand
                            if ~isempty(strfind(Words{k,4},'_p'))
                                EMDPP1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMDP1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end 
                        elseif ~isempty(strfind(Words{k,4},'_P'))
                            % Active Power
                            EMAP(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);                        
                        elseif ~isempty(strfind(Words{k,4},'_RE'))
                           % Reactive Energy  
                            if ~isempty(strfind(Words{k,4},'_parcial'))
                                EMREP1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMRE1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_E'))
                            % Active Energy
                            if ~isempty(strfind(Words{k,4},'_partial'))
                                EMAEP1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMAE1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end     
                        elseif ~isempty(strfind(Words{k,4},'_F'))
                            % Frequency
                            EMF1P(1:ndtEM1P(ind),ind)=eval([nomevar,'.',sensors{j}]);       
                        end
                    end
                    k=k+1;
                end
            end
        end
        ind=ind+1;
        eval('cd ..')
    end
end

disp('EM 340');

ndtEM=zeros(nEM,1);
dtEM=datetime(zeros(maxsamplesdt,nEM),'ConvertFrom','datenum');
EMTHDV_sys=zeros(maxsamplesdt,nEM);
EMTHDV_L1=zeros(maxsamplesdt,nEM);
EMTHDV_L2=zeros(maxsamplesdt,nEM);
EMTHDV_L3=zeros(maxsamplesdt,nEM);
EMVL_L_sys=zeros(maxsamplesdt,nEM);
EMVL_N_sys=zeros(maxsamplesdt,nEM);
EMVL1_L2=zeros(maxsamplesdt,nEM);
EMVL2_L3=zeros(maxsamplesdt,nEM);
EMVL3_L1=zeros(maxsamplesdt,nEM);
EMVL1_N=zeros(maxsamplesdt,nEM);
EMVL2_N=zeros(maxsamplesdt,nEM);
EMVL3_N=zeros(maxsamplesdt,nEM);
EMI_L1=zeros(maxsamplesdt,nEM);
EMI_L2=zeros(maxsamplesdt,nEM);
EMI_L3=zeros(maxsamplesdt,nEM);
EMTHDI_L1=zeros(maxsamplesdt,nEM);
EMTHDI_L2=zeros(maxsamplesdt,nEM);
EMTHDI_L3=zeros(maxsamplesdt,nEM);
EMPF_L1=zeros(maxsamplesdt,nEM);
EMPF_L2=zeros(maxsamplesdt,nEM);
EMPF_L3=zeros(maxsamplesdt,nEM);
EMPF_sys=zeros(maxsamplesdt,nEM);
EMApP_L1=zeros(maxsamplesdt,nEM);
EMApP_L2=zeros(maxsamplesdt,nEM);
EMApP_L3=zeros(maxsamplesdt,nEM);
EMApP_sys=zeros(maxsamplesdt,nEM);
EMRP_L1=zeros(maxsamplesdt,nEM);
EMRP_L2=zeros(maxsamplesdt,nEM);
EMRP_L3=zeros(maxsamplesdt,nEM);
EMRP_sys=zeros(maxsamplesdt,nEM);
EMAP_L1=zeros(maxsamplesdt,nEM);
EMAP_L2=zeros(maxsamplesdt,nEM);
EMAP_L3=zeros(maxsamplesdt,nEM);
EMAP_sys=zeros(maxsamplesdt,nEM);
EMRET=zeros(maxsamplesdt,nEM);
EMREP=zeros(maxsamplesdt,nEM);
EMDPP=zeros(maxsamplesdt,nEM);
EMDP=zeros(maxsamplesdt,nEM);
EMAET=zeros(maxsamplesdt,nEM);
EMAEP=zeros(maxsamplesdt,nEM);
EMAE_L1=zeros(maxsamplesdt,nEM);
EMAE_L2=zeros(maxsamplesdt,nEM);
EMAE_L3=zeros(maxsamplesdt,nEM);
EMF=zeros(maxsamplesdt,nEM);

for ind=1:nEM
    if IndEM(ind)~=0
        nome=Files{IndEM(ind)};
        eval(['cd ',nome]);
        d=dir;
        dev=d(3).name(1:2);
        nomef=[dev, '_',nome];
        eval(['load ',nomef]);
        nomevar=['data_',nomef];
        if nome(6)=='C' % asneira do Sergio
            nome(6)='c';
            nomef=[dev, '_',nome];
            nomevar=['data_',nomef];
            ind=2;
        else
            ind=1;
        end
        
        eval(['instant=',nomevar,'.instant;']);
        ndtEM(ind)=length(instant);
    
        dtEM(1:ndtEM(ind),ind)=datetime(instant);
        
        eval(['sensors=fieldnames(',nomevar,');']);
        nsen=length(sensors);
     
        for j=1:nsen
            if ~isempty(strmatch(sensors(j),'instant'))==0
                k=1;
                while k<=ntxt
                    if ~isempty(strmatch(sensorlist(k,:),sensors{j}))
                        if ~isempty(strfind(Words{k,4},'_V'))              
                            % voltages
                            if ~isempty(strfind(Words{k,4},'THD'))   
                                %Harmonic Distortion
                                if ~isempty(strfind(Words{k,4},'sys'))   
                                    EMTHDV_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                                elseif ~isempty(strfind(Words{k,4},'L1'))
                                    EMTHDV_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                                elseif ~isempty(strfind(Words{k,4},'L2'))
                                    EMTHDV_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                                elseif ~isempty(strfind(Words{k,4},'L3'))
                                    EMTHDV_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                                end
                            elseif ~isempty(strfind(Words{k,4},'sys'))
                                if ~isempty(strfind(Words{k,4},'L_L'))
                                    EMVL_L_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                                else
                                    EMVL_N_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);       
                                end
                            elseif ~isempty(strfind(Words{k,4},'L1_L2'))
                                EMVL1_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2_L3'))
                                EMVL2_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);                            
                            elseif ~isempty(strfind(Words{k,4},'L3_L1'))
                                EMVL3_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);      
                            elseif ~isempty(strfind(Words{k,4},'L1'))
                                EMVL1_N(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);                              
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMVL2_N(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);   
                            else
                                EMVL3_N(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);    
                            end
                            
                        elseif ~isempty(strfind(Words{k,4},'_I'))
                            % Current
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMI_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMI_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMI_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end
                        elseif ~isempty(strfind(Words{k,4},'_C'))
                            % Current THD
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMTHDI_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMTHDI_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMTHDI_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end     
                        elseif ~isempty(strfind(Words{k,4},'_PF'))
                            % Power Factor
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMPF_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMPF_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                             elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMPF_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            else
                                EMPF_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_ApP'))
                            % Apparent Power
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMApP_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMApP_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                             elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMApP_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            else
                                EMApP_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_RP'))
                            % Reactive Power
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMRP_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMRP_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                             elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMRP_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            else
                                EMRP_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_P'))
                            % Active Power
                            if ~isempty(strfind(Words{k,4},'L1'))
                                EMAP_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMAP_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                             elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMAP_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            else
                                EMAP_sys(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end                        
                        elseif ~isempty(strfind(Words{k,4},'_RE'))
                           % Reactive Energy  
                            if ~isempty(strfind(Words{k,4},'_total'))
                                EMRET(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMREP(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_ED'))
                           % Energy Demand
                            if ~isempty(strfind(Words{k,4},'_p'))
                                EMDPP(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                EMDP(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end 
                        elseif ~isempty(strfind(Words{k,4},'_E'))
                            % Active Energy
                            if ~isempty(strfind(Words{k,4},'_total'))
                                EMAET(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'_partial'))
                                EMAEP(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L1'))
                                EMAE_L1(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                EMAE_L2(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                             elseif ~isempty(strfind(Words{k,4},'L3'))
                                EMAE_L3(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end     
                        elseif ~isempty(strfind(Words{k,4},'_F'))
                           % Reactive Energy                       
                           EMF(1:ndtEM(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        end
                    end
                    k=k+1;
                end
            end
        end
        ind=ind+1;
        eval('cd ..')
    end
end

disp('Inverter');

maxsamplesINV=24*60;
ndtINV=zeros(nINV,1);
dtINV=datetime(zeros(maxsamplesINV,nINV),'ConvertFrom','datenum');
INVGE=zeros(maxsamplesINV,nINV);
INVGP=zeros(maxsamplesINV,nINV);
INVBCDP=zeros(maxsamplesINV,nINV);
INVPMI_L1=zeros(maxsamplesINV,nINV);
INVPMI_L2=zeros(maxsamplesINV,nINV);
INVPMI_L3=zeros(maxsamplesINV,nINV);
INVPMAP_sys=zeros(maxsamplesINV,nINV);
INVPMAP_L1=zeros(maxsamplesINV,nINV);
INVPMAP_L2=zeros(maxsamplesINV,nINV);
INVPMAP_L3=zeros(maxsamplesINV,nINV);
INVPMRP_sys=zeros(maxsamplesINV,nINV);
INVPMRP_L1=zeros(maxsamplesINV,nINV);
INVPMRP_L2=zeros(maxsamplesINV,nINV);
INVPMRP_L3=zeros(maxsamplesINV,nINV);
INVPMApP_sys=zeros(maxsamplesINV,nINV);
INVPMApP_L1=zeros(maxsamplesINV,nINV);
INVPMApP_L2=zeros(maxsamplesINV,nINV);
INVPMApP_L3=zeros(maxsamplesINV,nINV);
INVPMV_L1=zeros(maxsamplesINV,nINV);
INVPMV_L2=zeros(maxsamplesINV,nINV);
INVPMV_L3=zeros(maxsamplesINV,nINV);
INVPMPF=zeros(maxsamplesINV,nINV);
INVPMF=zeros(maxsamplesINV,nINV);
INVDCV_L1=zeros(maxsamplesINV,nINV);
INVDCV_L2=zeros(maxsamplesINV,nINV);
INVDCV_L3=zeros(maxsamplesINV,nINV);
INVDCI_L1=zeros(maxsamplesINV,nINV);
INVDCI_L2=zeros(maxsamplesINV,nINV);
INVDCI_L3=zeros(maxsamplesINV,nINV);
INVDCP_L1=zeros(maxsamplesINV,nINV);
INVDCP_L2=zeros(maxsamplesINV,nINV);
INVDCP_L3=zeros(maxsamplesINV,nINV);
INVDCP_sys=zeros(maxsamplesINV,nINV);
INVECB=zeros(maxsamplesINV,nINV);
INVECG=zeros(maxsamplesINV,nINV);
INVECPV=zeros(maxsamplesINV,nINV);
INVCR=zeros(maxsamplesINV,nINV);
INVECT=zeros(maxsamplesINV,nINV);
INVPCG=zeros(maxsamplesINV,nINV);
INVPCB=zeros(maxsamplesINV,nINV);
INVPCPV=zeros(maxsamplesINV,nINV);
INVYD=zeros(maxsamplesINV,nINV);
INVYM=zeros(maxsamplesINV,nINV);
INVYY=zeros(maxsamplesINV,nINV);
INVYT=zeros(maxsamplesINV,nINV);
INVMS=zeros(maxsamplesINV,nINV);
INVIS=zeros(maxsamplesINV,nINV);
INVPL=zeros(maxsamplesINV,nINV);
INVWT=zeros(maxsamplesINV,nINV);
INVBSC=zeros(maxsamplesINV,nINV);
INVBCC=zeros(maxsamplesINV,nINV);
INVBDC=zeros(maxsamplesINV,nINV);
INVBGC=zeros(maxsamplesINV,nINV);
INVBT=zeros(maxsamplesINV,nINV);
INVBNC=zeros(maxsamplesINV,nINV); 
INVBV=zeros(maxsamplesINV,nINV); 
INVBASC=zeros(maxsamplesINV,nINV); 
INVI_L1=zeros(maxsamplesINV,nINV); 
INVI_L2=zeros(maxsamplesINV,nINV); 
INVI_L3=zeros(maxsamplesINV,nINV); 
INVAP_sys=zeros(maxsamplesINV,nINV);
INVAP_L1=zeros(maxsamplesINV,nINV);
INVAP_L2=zeros(maxsamplesINV,nINV);
INVAP_L3=zeros(maxsamplesINV,nINV);
INVRP_sys=zeros(maxsamplesINV,nINV);
INVApP_sys=zeros(maxsamplesINV,nINV);
INVV_L1=zeros(maxsamplesINV,nINV);
INVV_L2=zeros(maxsamplesINV,nINV);
INVV_L3=zeros(maxsamplesINV,nINV);
INVPF=zeros(maxsamplesINV,nINV);
INVF=zeros(maxsamplesINV,nINV);

INVBatAPSP=zeros(maxsamplesINV,nINV);
INVBatRPSP=zeros(maxsamplesINV,nINV);
INVDeltaCos=zeros(maxsamplesINV,nINV);
INVTotDCcharge=zeros(maxsamplesINV,nINV);
INVTotDCdischarge=zeros(maxsamplesINV,nINV);
INVTotACcharge=zeros(maxsamplesINV,nINV);
INVTotACdischarge=zeros(maxsamplesINV,nINV);
INVTotDCPV=zeros(maxsamplesINV,nINV);
INVTotDCPV1=zeros(maxsamplesINV,nINV);
INVTotDCPV2=zeros(maxsamplesINV,nINV);
INVTotDCPV3=zeros(maxsamplesINV,nINV);
INVTotACenergy=zeros(maxsamplesINV,nINV);
INVTotDCpower=zeros(maxsamplesINV,nINV);
INVTotACchargegrid=zeros(maxsamplesINV,nINV);

for ind=1:nINV
    nome=Files{IndINV(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    new=~isempty(strfind(d(3).name,'new'));
    newnew=~isempty(strfind(d(3).name,'newvars'));
    control=~isempty(strfind(d(3).name,'control'));
    if control
        nomef=[dev, '_plenticore_kostal_nv_control'];
    elseif newnew
        nomef=[dev, '_plenticore_kostal_newvars'];
    elseif new
        nomef=[dev, '_plenticore_kostal_new'];
    else
        nomef=[dev, '_plenticore_kostal'];
    end
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    if control
        newnew=1;
    end
    if newnew
        new=1;
    end
        
    eval(['instant=',nomevar,'.instant;']);
    ndtINV(ind)=length(instant);

    dtINV(1:ndtINV(ind),ind)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'_inv')) 
                        if ~isempty(strfind(Words{k,4},'_e'))
                            INVGE(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        else
                            INVGP(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        end
                    end
                    if ~isempty(strfind(Words{k,4},'bat_')) 
                        INVBCDP(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    end
                    if ~isempty(strfind(Words{k,4},'pc_pm'))              
                        % powermeter
                        if ~isempty(strfind(Words{k,4},'_I'))   
                            %Current
                            if ~isempty(strfind(Words{k,4},'L1'))
                                INVPMI_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVPMI_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVPMI_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            end
                        elseif ~isempty(strfind(Words{k,4},'_AP'))
                            %Active Power
                            if ~isempty(strfind(Words{k,4},'T_AP'))
                                INVPMAP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L1'))
                                INVPMAP_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);      
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVPMAP_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVPMAP_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end
                        elseif ~isempty(strfind(Words{k,4},'_RP'))
                            %Reactive Power
                            if ~isempty(strfind(Words{k,4},'T_RP'))
                                INVPMRP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L1'))
                                INVPMRP_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);      
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVPMRP_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVPMRP_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end
                        elseif ~isempty(strfind(Words{k,4},'_ApP'))
                            %Apparent Power
                            if ~isempty(strfind(Words{k,4},'T_ApP'))
                                INVPMApP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L1'))
                                INVPMApP_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);      
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVPMApP_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVPMApP_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end                            
                        elseif ~isempty(strfind(Words{k,4},'_V'))
                            % Voltage
                            if ~isempty(strfind(Words{k,4},'L1'))
                                INVPMV_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVPMV_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVPMV_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end
                        elseif ~isempty(strfind(Words{k,4},'_cos'))
                            % Power Factor
                            INVPMPF(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_freq'))
                            % Frequency
                            INVPMF(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);        
                        end        

                    elseif ~isempty(strfind(Words{k,4},'_DC'))
                        % DC Variables
                        if ~isempty(strfind(Words{k,4},'_V'))
                            % Voltage
                            if ~isempty(strfind(Words{k,4},'L1'))
                                INVDCV_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVDCV_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                INVDCV_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end     
                        elseif ~isempty(strfind(Words{k,4},'_I'))
                            % Current
                            if ~isempty(strfind(Words{k,4},'L1'))
                                INVDCI_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVDCI_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            else
                                INVDCI_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end     
                        elseif ~isempty(strfind(Words{k,4},'_P'))
                            % Power
                            if ~isempty(strfind(Words{k,4},'L1'))
                                INVDCP_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L2'))
                                INVDCP_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                            elseif ~isempty(strfind(Words{k,4},'L3'))
                                INVDCP_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            elseif ~isempty(strfind(Words{k,4},'_T_'))
                                INVDCP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                            end 
                        end
                    %Others

                    elseif ~isempty(strfind(Words{k,4},'_TH')) 
                        % Home Energy Consumption
                        if ~isempty(strfind(Words{k,4},'_Bat'))
                            % from Battery
                            INVECB(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_Gri'))
                            % from Grid
                            INVECG(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_pv'))
                            % from PV
                            INVECPV(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_C_r'))
                            % Total Consumption Rate
                            INVCR(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        else
                            INVECT(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        end
                    elseif ~isempty(strfind(Words{k,4},'_H_P_')) 
                        % Home Power Consumption
                        if ~isempty(strfind(Words{k,4},'_Bat'))
                            % from Battery
                            INVPCB(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_Gri'))
                            % from Grid
                            INVPCG(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_pv'))
                            % from PV
                            INVPCPV(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        end
                    elseif ~isempty(strfind(Words{k,4},'_E_Y')) 
                        % Yield
                        if ~isempty(strfind(Words{k,4},'_D'))
                            % daily
                            INVYD(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_M'))
                            % monthly
                            INVYM(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_Y_E_Y'))
                            % yearly
                            INVYY(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                       elseif ~isempty(strfind(Words{k,4},'_T'))
                            % yearly
                            INVYT(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        end
                    elseif ~isempty(strfind(Words{k,4},'_E_M_S'))
                        INVMS(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_Inv_S'))
                        INVIS(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_P_L_'))
                        INVPL(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_wt'))
                        INVWT(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_Bat'))   
                        % Battery
                        if ~isempty(strfind(Words{k,4},'_S_C'))
                            INVBSC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_C_I'))
                            INVBCC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'_CD_I'))
                            INVBDC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'Cap'))
                            INVBGC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);  
                        elseif ~isempty(strfind(Words{k,4},'_temp'))
                            INVBT(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        elseif ~isempty(strfind(Words{k,4},'_N_C'))
                            INVBNC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        elseif ~isempty(strfind(Words{k,4},'_V'))
                            INVBV(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        elseif ~isempty(strfind(Words{k,4},'_SOC'))
                            INVBASC(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        end
                    elseif ~isempty(strfind(Words{k,4},'_I'))   
                        %Current
                        if ~isempty(strfind(Words{k,4},'L1'))
                            INVI_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            INVI_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            INVI_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);                           
                        end
                    
                    elseif ~isempty(strfind(Words{k,4},'_AP'))
                        %Active Power
                        if ~isempty(strfind(Words{k,4},'T_AP'))
                            INVAP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        elseif ~isempty(strfind(Words{k,4},'L1'))
                            INVAP_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);      
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            INVAP_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            INVAP_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        end
                    elseif ~isempty(strfind(Words{k,4},'_RP'))
                        %Reactive Power
                        INVRP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_ApP'))
                        %Apparent Power
                        INVApP_sys(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);                     
                    elseif ~isempty(strfind(Words{k,4},'_V'))
                        % Voltage
                        if ~isempty(strfind(Words{k,4},'L1'))
                            INVV_L1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L2'))
                            INVV_L2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'L3'))
                            INVV_L3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]); 
                        end
                    elseif ~isempty(strfind(Words{k,4},'_cos'))
                        % Power Factor
                        INVPF(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'_freq'))
                        % Frequency
                        INVF(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);        
                    elseif ~isempty(strfind(Words{k,4},'pc_ap_setp'))
                        % Active Power Setpoint
                        INVBatAPSP(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_rp_setp'))
                        % Reactive Power Setpoint
                        INVBatRPSP(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_dcosphi_setp'))
                        % Delta Cos Phi SetPoint
                        INVDeltaCos(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_ce'))
                        % Total DC Charge Energy (DC-side-to-battery)
                        INVTotDCcharge(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_de'))
                        % Total DC Charge Energy (DC-side-to-battery)
                        INVTotDCdischarge(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tac_ce'))
                        % Total AC Charge Energy (AC-side-to-battery)
                        INVTotACcharge(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tac_de'))
                        % Total AC Discharge Energy (AC-side-from-battery)
                        INVTotACdischarge(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_pve'))
                        % Total DC energy (Sum of all inputs)
                        INVTotDCPV(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_epv1'))
                        % Total DC energy from PV1
                        INVTotDCPV1(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_epv2'))
                        % Total DC energy from PV2
                        INVTotDCPV2(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_epv3'))
                        % Total DC energy from PV3
                        INVTotDCPV3(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_teac_grid'))
                        % Total AC energy AC_side to grid
                        INVTotACenergy(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tdc_p'))
                        % Total DC power
                        INVTotDCpower(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'pc_tac_ceg'))
                        % Total AC charge energy grid
                        INVTotACchargegrid(1:ndtINV(ind),ind)=eval([nomevar,'.',sensors{j}]);
                    end
                    
                end
                k=k+1;
            end
        end
    end
    ind=ind+1;
    eval('cd ..')
end

disp('Smart Plugs');
ndtSP=zeros(nSP,1);
dtSP=datetime(zeros(maxsamplesdt,nSP),'ConvertFrom','datenum');
SPV=zeros(maxsamplesdt,nSP);
SPI=zeros(maxsamplesdt,nSP);
SPAP=zeros(maxsamplesdt,nSP);
SPAE=zeros(maxsamplesdt,nSP);
SPRssi=zeros(maxsamplesdt,nSP);
SPOn=zeros(maxsamplesdt,nSP);

for ind=1:nSP
    if TabSP(ind)
        i=ind;
        nome=Files{IndSP(ind)};
        eval(['cd ',nome]);
        d=dir;
        dev=d(3).name(1:2);
        nomef=[dev, '_',nome];
        eval(['load ',nomef]);
        nomevar=['data_',nomef];
        
        eval(['instant=',nomevar,'.instant;']);
        [ndtSP(i),~]=size(instant);
        dtSP(1:ndtSP(i),i)=datetime(instant);
        
        eval(['sensors=fieldnames(',nomevar,');']);
        nsen=length(sensors);
        for j=1:nsen
            if length(strmatch(sensors(j),'instant'))==0
                k=1;
                while k<=ntxt
                    if length(strmatch(sensorlist(k,:),sensors{j}))
                        if ~isempty(strfind(Words{k,4},'Voltage'))
                            SPV(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'Current'))
                            SPI(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'Power'))
                            SPAP(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'Total'))
                            SPAE(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);
                        elseif ~isempty(strfind(Words{k,4},'Rssi'))
                            SPRssi(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);    
                        elseif ~isempty(strfind(Words{k,4},'On'))
                            SPOn(1:ndtSP(i),i)=eval([nomevar,'.',sensors{j}]);
                        end
                    end
                    k=k+1;
                end
            end
        end
        eval('cd ..')
    end
end

disp('Weather Station');
ndtWS=zeros(nWS,1);
dtWS=datetime(zeros(maxsamplesINV,nWS),'ConvertFrom','datenum');
WS_RAD=zeros(maxsamplesINV,nWS);
WS_AT=zeros(maxsamplesINV,nWS);
WS_RH=zeros(maxsamplesINV,nWS);


for ind=1:nWS
    i=ind;
    nome=Files{IndWS(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    ndtWS(i)=length(instant);
    dtWS(1:ndtWS(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'ws_rad'))
                        WS_RAD(1:ndtWS(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ws_temp'))
                        WS_AT(1:ndtWS(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ws_hum'))
                        WS_RH(1:ndtWS(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end


disp('Hall 1');
maxsamplesSPWS=24*60*2;
ndtH1=zeros(nSPWS_hall1,1);
dtH1=datetime(zeros(maxsamplesSPWS,nSPWS_hall1),'ConvertFrom','datenum');
H1_AT=zeros(maxsamplesSPWS,nSPWS_hall1);
H1_RH=zeros(maxsamplesSPWS,nSPWS_hall1);
H1_L=zeros(maxsamplesSPWS,nSPWS_hall1);

for ind=1:nSPWS_hall1
    i=ind;
    nome=Files{IndSPWS_hall1(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    [ndtH1(i), ~]=size(instant);
    dtH1(1:ndtH1(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'em_temp_c'))
                        H1_AT(1:ndtH1(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_hum_c'))
                        H1_RH(1:ndtH1(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_light_c'))
                        H1_L(1:ndtH1(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

disp('bedroom_1_2');
ndtB12=zeros(nSPWS_bedroom_1_2,1);
dtB12=datetime(zeros(maxsamplesSPWS,nSPWS_bedroom_1_2),'ConvertFrom','datenum');
B12_AT=zeros(maxsamplesSPWS,nSPWS_bedroom_1_2);
B12_RH=zeros(maxsamplesSPWS,nSPWS_bedroom_1_2);
B12_L=zeros(maxsamplesSPWS,nSPWS_bedroom_1_2);

for ind=1:nSPWS_bedroom_1_2
    i=ind;
    nome=Files{IndSPWS_bedroom_1_2(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    ndtB12(i)=size(instant,1);
    dtB12(1:ndtB12(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'em_temp_qp'))
                        B12_AT(1:ndtB12(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_hum_qp'))
                        B12_RH(1:ndtB12(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_light_qp'))
                        B12_L(1:ndtB12(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

disp('bedroom_1_4');
ndtB14=zeros(nSPWS_bedroom_1_4,1);
dtB14=datetime(zeros(maxsamplesSPWS,nSPWS_bedroom_1_4),'ConvertFrom','datenum');
B14_AT=zeros(maxsamplesSPWS,nSPWS_bedroom_1_4);
B14_RH=zeros(maxsamplesSPWS,nSPWS_bedroom_1_4);
B14_M=zeros(maxsamplesSPWS,nSPWS_bedroom_1_4);
B14_WT=zeros(maxsamplesSPWS,nSPWS_bedroom_1_4);

for ind=1:nSPWS_bedroom_1_4
    i=ind;
    nome=Files{IndSPWS_bedroom_1_4(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    [ndtB14(i),~]=size(instant);
    dtB14(1:ndtB14(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'em_temp_qg'))
                        B14_AT(1:ndtB14(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_hum_qg'))
                        B14_RH(1:ndtB14(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_mov_qg'))
                        B14_M(1:ndtB14(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_wt_qg'))
                        B14_WT(1:ndtB14(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

disp('Lounge');
ndtL=zeros(nSPWS_lounge,1);
dtL=datetime(zeros(maxsamplesSPWS,nSPWS_lounge),'ConvertFrom','datenum');
L_AT=zeros(maxsamplesSPWS,nSPWS_lounge);
L_RH=zeros(maxsamplesSPWS,nSPWS_lounge);
L_M=zeros(maxsamplesSPWS,nSPWS_lounge);

for ind=1:nSPWS_lounge
    i=ind;
    nome=Files{IndSPWS_lounge(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    [ndtL(i),~]=size(instant);
    dtL(1:ndtL(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'em_temp_s'))
                        L_AT(1:ndtL(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_hum_s'))
                        L_RH(1:ndtL(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'em_mov_s'))
                        L_M(1:ndtL(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

disp('AC');
ndtAC=zeros(nAC,1);
dtAC=datetime(zeros(maxsamplesINV,nAC),'ConvertFrom','datenum');
AC_PS=zeros(maxsamplesINV,nAC);
AC_SM=zeros(maxsamplesINV,nAC);
AC_EM=zeros(maxsamplesINV,nAC);
AC_TM=zeros(maxsamplesINV,nAC);
AC_OM=zeros(maxsamplesINV,nAC);
AC_RT=zeros(maxsamplesINV,nAC);
AC_IT=zeros(maxsamplesINV,nAC);
AC_OT=zeros(maxsamplesINV,nAC);
AC_FS=zeros(maxsamplesINV,nAC);


for ind=1:nAC
    i=ind;
    nome=Files{IndAC(ind)};
    eval(['cd ',nome]);
    d=dir;
    dev=d(3).name(1:2);
    nomef=[dev, '_',nome];
    eval(['load ',nomef]);
    nomevar=['data_',nomef];
    
    eval(['instant=',nomevar,'.instant;']);
    [ndtAC(i),~]=size(instant);
    dtAC(1:ndtAC(i),i)=datetime(instant);
    
    eval(['sensors=fieldnames(',nomevar,');']);
    nsen=length(sensors);
    for j=1:nsen
        if length(strmatch(sensors(j),'instant'))==0
            k=1;
            while k<=ntxt
                if length(strmatch(sensorlist(k,:),sensors{j}))
                    if ~isempty(strfind(Words{k,4},'ac_ps'))
                        AC_PS(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ac_sm'))
                        AC_SM(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ac_em'))
                        AC_EM(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ac_tm'))
                        AC_TM(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);                        
                    elseif ~isempty(strfind(Words{k,4},'ac_opm'))
                        AC_OM(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);                        
                    elseif ~isempty(strfind(Words{k,4},'ac_setp'))
                        AC_RT(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]); 
                    elseif ~isempty(strfind(Words{k,4},'ac_intemp'))
                        AC_IT(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ac_outemp'))
                        AC_OT(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    elseif ~isempty(strfind(Words{k,4},'ac_fsp'))
                        AC_FS(1:ndtAC(i),i)=eval([nomevar,'.',sensors{j}]);
                    end
                end
                k=k+1;
            end
        end
    end
    eval('cd ..')
end

save(dir_name,'ndt','dt','V','I','F','AP','RP','ApP','PF','AE','IRE','CRE','ndt3','dt3','average_F','T_AP','T_RP','T_ApP','T_PF','T_AE','T_IRE','T_CRE');
if nEM>0
    save(dir_name,'EMVL1_L2', 'EMVL1_N','EMVL2_L3', 'EMVL2_N', 'EMVL3_L1', 'EMVL3_N', 'EMVL_N_sys', 'EMVL_L_sys', 'EMI_L1', 'EMI_L3', 'EMAP_sys', 'EMAP_L1', 'EMAP_L2', ...
    'EMAP_L3', 'EMApP_sys', 'EMApP_L1', 'EMApP_L2', 'EMApP_L3', 'EMRP_sys', 'EMRP_L1', 'EMRP_L2', 'EMRP_L3', 'EMPF_sys', 'EMPF_L1', 'EMPF_L2', 'EMPF_L3', ...
    'EMTHDV_sys', 'EMTHDV_L1', 'EMTHDV_L2', 'EMTHDV_L3', 'EMTHDI_L1', 'EMTHDI_L2', 'EMTHDI_L3', 'EMF', 'EMAE_L1', 'EMAE_L2', 'EMAE_L3', 'EMAET', 'EMAEP', 'EMRET', 'EMREP', ...
    'EMDP', 'EMDPP', 'dtEM', 'ndtEM', '-append');
    if exist('EMI_L2','var')
        save(dir_name,'EMI_L2','-append')
    end
end
if nEM1P>0
    save(dir_name,'dtEM1P','ndtEM1P','EMV','EMI','EMAP','EMApP','EMRP','EMPF','EMF1P','EMAE1P','EMAEP1P','EMRE1P','EMREP1P','EMDP1P','EMDPP1P','dtEM1P', 'ndtEM1P','-append');
end

if nINV>0
    save(dir_name,'INVPMPF','INVPMF','INVPMI_L1','INVPMI_L2','INVPMI_L3','INVPMAP_sys','INVPMAP_L1','INVPMAP_L2','INVPMAP_L3','INVPMRP_sys','INVPMRP_L1','INVPMRP_L2','INVPMRP_L3','INVPMApP_sys','INVPMApP_L1','INVPMApP_L2','INVPMApP_L3',...
        'INVPMV_L1','INVPMV_L2','INVPMV_L3','INVDCI_L1','INVDCI_L2','INVDCI_L3','INVDCV_L1', 'INVDCV_L2','INVDCV_L3','INVDCP_L1','INVDCP_L2','INVDCP_L3','INVPF', 'INVF', 'INVI_L1', 'INVI_L2', 'INVI_L3', 'INVAP_sys', 'INVAP_L1', ...
        'INVAP_L2', 'INVAP_L3', 'INVRP_sys', 'INVApP_sys', 'INVV_L1', 'INVV_L2', 'INVV_L3', 'INVBSC','INVBCC', 'INVBDC', 'dtINV', 'ndtINV', '-append');
    if new
        save(dir_name, 'INVDCP_sys', 'INVBGC', 'INVBT', 'INVBNC', 'INVBV', 'INVECB', 'INVECG', 'INVECPV', 'INVECT', 'INVPCB', 'INVPCG', 'INVPCPV', 'INVYD', 'INVYM', 'INVYY', 'INVYT', 'INVMS', 'INVIS', 'INVPL', 'INVWT', 'INVCR','-append')
        if exist('INVBASC','var')
            save(dir_name,'INVBASC','-append')
        end
        if newnew
            save(dir_name, 'INVGP','INVGE', 'INVBCDP', '-append')
            if control
                save(dir_name, 'INVBatAPSP','INVBatRPSP','INVDeltaCos','INVTotDCcharge','INVTotDCdischarge','INVTotACcharge','INVTotACdischarge','INVTotDCPV','INVTotDCPV1','INVTotDCPV2','INVTotDCPV3','INVTotACenergy','INVTotDCpower','INVTotACchargegrid','-append')
            end
        end
    end
end

if nSP>0
    save(dir_name,'SPV', 'SPI', 'SPAP', 'SPAE', 'SPRssi', 'SPOn', 'dtSP', 'ndtSP', '-append');
end

if nWS>0
    save(dir_name,'WS_RAD', 'WS_AT', 'WS_RH', 'dtWS', 'ndtWS', '-append');
end
 
if nSPWS_hall1>0
    save(dir_name,'H1_L', 'H1_AT', 'H1_RH', 'dtH1', 'ndtH1', '-append');
end

if nSPWS_bedroom_1_2>0
    save(dir_name,'B12_L', 'B12_AT', 'B12_RH', 'dtB12', 'ndtB12', '-append');
end

if nSPWS_bedroom_1_4>0
    save(dir_name,'B14_M', 'B14_WT','B14_AT', 'B14_RH', 'dtB14', 'ndtB14', '-append');
end

if nSPWS_lounge>0
    save(dir_name,'L_M','L_AT', 'L_RH', 'dtL', 'ndtL', '-append');
end

if nAC>0
    save(dir_name,'AC_PS','AC_SM', 'AC_EM', 'AC_TM', 'AC_OM', 'AC_RT', 'AC_IT', 'AC_OT', 'AC_FS', 'dtAC', 'ndtAC', '-append');
end


eval('cd ..')

