function result = STHEmass(HXdata,Fixed,outputfname,TEMA,plotti)  %(do,tb_arr,Lpass,Ntb_pass,DOTL,Nbfl,stb,P)
% HXdata
% Fixed
% outputfname
% TEMA: character array of 3:
%       TEMA(1): front-end head
%       TEMA(2): shell
%       TEMA(3): rear-end head
% plotti
%
% SHELL THICKNESS: 
% 1) if only shell inner D given, sizing is performed against internal
% pressure and elastic&plastic buckling and internal pressure 
% --> s_sh, D_shO obtained.
% 2) if HXdata contains shell inner&outer D's and thickness, this is
% compared against sizing results & TEMA-specified minimum shell thickness 
% --> error given is inadequate.

fHead=TEMA(1); 
Shell=TEMA(2);
rHead=TEMA(3);

metalprops = Fixed.metalprops;

fubar=0;                                    % default: everything ok to continue.
switch fHead
    case '0'                                % no front head: for DH SMR
        Nflng = 0;                          % Nflng = number of flanges
    case 'A'
        Nflng = 3;
    case 'B'
        Nflng = 2;
    case 'C'
        Nflng = 3;  
    case 'N'
        Nflng = 1; 
    otherwise
        fubar=1;
end
RearHeadSmallFlange=0;

switch rHead
    case '0'                                % no front head: for DH SMR
        Nflng = Nflng+0;
    case 'L'
        Nflng = Nflng+3;
    case 'M'
        Nflng = Nflng+1; 
    case 'N'  
        Nflng = Nflng+1;  
    case 'T'   
        Nflng = Nflng+2;   
        RearHeadSmallFlange = RearHeadSmallFlange+1;
    case 'U'   
        Nflng = Nflng+0;
    otherwise
        fubar=1;
end

s_TbPassDivPlate = 0.016;                   % tube pass division plate in the shell: 16 mm [HEDH 4.2.5-3]

TbPasses = HXdata.Npasses;
do=HXdata.do_mm;
stb=HXdata.s_mm;
P=HXdata.P_mm;
D_OTL = HXdata.D_OTL;                       % [m] bundle diameter, outer tube limit [m]         
D_ShIn=HXdata.D_ShIn;                       % [m] shell inside diameter  
LStraight = HXdata.LStraight;               % [m] Lstraight = straight tube in pass INCLUDING covered by baffles, NOT INCLUDING inside tubesheet
NtbPerPass = HXdata.NtbPerPass;             % number of tubes per pass
N_bfls=HXdata.N_bfls;
N_HolesPerBaffle = HXdata.N_HolesPerBaffle;
AbflCond=HXdata.A_BflFr;                       
SStripWidth =HXdata.SStripWidth;
Lshell = LStraight;                         % [m] shell length = straight tube length     

%ptbmin=Fixed.sizingdata.ptbmin;            % tube-side min p [MPa]
ptbmax=Fixed.sizingdata.ptbmax;             % tube-side max p [MPa]
pshmin=Fixed.sizingdata.pshmin;           	% shell-side min p [MPa]
pshmax=Fixed.sizingdata.pshmax;          	% shell-side max p [MPa]                            % Ttbmax=fscanf(fid,'%g', [1]);     

if (pshmin<0)                              	% p_out_s = max overpressure on the shell from outside [MPa] (always >0)
    p_out_s = abs(pshmin);                  
else                        
    p_out_s = 0;
end                                                                                            %ptblomm(1) = pshmax-ptbmin; %ptblomm(2) = pshmax; %ptblomm(3) = .1-ptbmin; %ptblomm(4) = 0;                             % p_out_tb = max(ptblomm);                    % p_out_tb = max overpressure on the shell from outside [MPa] (always >=0)% Lshell = 0.001*Lshell;                    % shell length [m] <- [mm]s_avtb = 0.004;

s_tb = 0.001*stb;                           % tube wall thickness [m] <- [mm]
do = 0.001*do;                              % tube outer diameter  [m] <- [mm]
di = do - (2*s_tb);
if (isfield(metalprops,'rho'))
    rho=metalprops.rho;                     % material [kg/m3]
else
    rho = 7850;                            	% material [kg/m3]
end
s_bfl = HXdata.s_bfl;                       % baffle plate thickness [m]  0.001*12.5;
Nrods = 8;                                  % number of tie rods
dRod = 0.001* 12.5;                         % diameter of rod 12.5 mm
N_tbplates = 2;                             % not U-tube so 2 tubesheets

% CHANNEL LENGTHS, TUBE PASSES, TUBE PLATE NUMBERS ACCORDING TO FIXED HEAD DESIGN
% ================================================================================
%
% LENGTHS AND THICKNESSES
% -----------------------
%
% Channel length and diameter
if ((D_OTL+0.050) < (D_ShIn/1.25))
    D_channel = D_OTL+0.050;
else
    D_channel = D_ShIn;
end
Atb=NtbPerPass*0.25*pi*(di^2);                                  % free-flow area into tubes in one pass [m^2]
Anzl=Atb/1.5;                                                   % assume 50% greater velocity ok for DH pipe in/out
Dnzl=(Anzl*4/pi)^.5;                                            % cold-side nozzle diametre                            % entry (and exit) head chamber length [m]
Lhead=max((2*Dnzl),(0.75*D_channel));                           % preliminary estimate for determining channel thickness

% thickness, outer D of channel
if ((fHead~='0') && (rHead~='0'))                               % channel sizing only if there is a channel
    s_channel = shellthickness(p_out_s,pshmax,(D_ShIn*1000),Lhead,0.5,50, metalprops); % s_shell [mm]
    s_channel = s_channel*0.001;                                % s_shell [m]
    Du_channel = D_channel + (2*s_channel);                     % shell O.D. [m]
else                                                            % if there's no channel, dimensions are zero.
    s_channel = 0;
    Du_channel = 0;
end

% thickness, outer D of shell.
s_shell = shellthickness(p_out_s,pshmax,(D_ShIn*1000),(Lshell*1000),0.5,50, metalprops); % s_shell sized against internal pressure and elast/plast buckling [mm]
s_shell_min = TEMA_minShellThickness(D_ShIn,2);                 % TEMA_minShellThickness(Dshell,material) -- material 0: C-steel/pipe, 1: C-steel/plate, 2: alloy 
s_shell = 0.001*max(s_shell,s_shell_min);                       % s_shell [m]: greater of TEMA minimum & pressure vessel sizing results 
if (isfield(HXdata,'s_sh'))                                     %
    if s_shell > HXdata.s_sh                                    % did TEMA minimum or pressure vessel code required a thicker-than-specified shell?
        fubar=1;                                                 % if yes - set error
        fprintf(1, '\n\n  INSUFFICIENT SHELL THICKNESS ');      % and print an error
        fprintf(1, 'SPECIFIED -- INCREASE PIPE SCHEDULE OR ');  % message
        fprintf(1, 'USE THICKER PLATE !!!');                                     
        fprintf(1, ' D = %5.3f m, L = %5.3f m, smin = %4.2f mm \n\n',D_ShIn,Lshell,s_shell*1000);
    else                                                        % otherwise everything a.o.k., shell specified at least as thick as min.required,
        s_shell = HXdata.s_sh;                                  % copy specified thickness to calculation thickness
        if (isfield(HXdata,'D_ShOut'))                          % check if shell outer D specified also, and if exists,
            Du_shell = HXdata.D_ShOut;                          % take the specified for calculation, after checking if 
            if abs(D_ShIn+2*s_shell-Du_shell)>0.0001            % if error between specified, calculated D_shO > 0.1 mm, 
                fubar=1;                                         % set error code (this should never happpen, but....)
            end
        else                                                    % If shell O.D. not specified, calculate the 
            Du_shell = D_ShIn + (2*s_shell);                    % shell O.D. [m] with shell I.D. and thickness
        end
    end
else                                                            % if no s_shell pre-specified, use the greater of TEMA minimum & pressure vessel sizing result
    Du_shell = D_ShIn + (2*s_shell);                            % with shell I.D. to calculate the shell O.D. [m] 
end

% thicknesses of tubesheet, channel cover, channel shell
s_plates = platethickness(ptbmax,pshmax,(D_ShIn*1000),(D_channel*1000),(s_shell*1000),(s_channel*1000),(do*1000),P,metalprops);
s_tbplate = s_plates(1)*0.001;                                  % tube plate thickness [m]
s_channelcover = s_plates(2)*0.001;                             % end cover thickness [m]
s_channelshell_int_p = 0.001*channelshell_intp(ptbmax,(D_ShIn*1000),metalprops);% s_internalp returns thickness against internal p [mm] (SFS-2611)
s_channelshell = max(s_shell,s_channelshell_int_p);             % channel shell thickness [m];
Du_channel = Du_channel + (2*s_channelshell);                 	% channel O.D. [m]

% flange, channel cover and tubesheet D 

if fubar==0
    Fdims = flangedimsSimple(max((Du_shell*1000),(Du_channel*1000)),(s_channel*1000),(s_tbplate*1000));
    D_Flng=Fdims.D;
    switch fHead
        case '0'
            D_fCoverplate=0;                        % no cover plate
            s_fChannelcover = 0;                    % no channel
            L_fhead=0;                              % no entry head chamber 
        case 'A'
            D_fCoverplate=D_Flng;
            s_fChannelcover = s_channelcover;       % flat plate bolted to flanges
            L_fhead=max((2*Dnzl),(0.75*D_channel)); % entry head chamber length
        case 'B'
            D_fCoverplate = Du_channel; 
            s_fChannelcover = s_channelshell;       % bonnet cover, estimate same thickness as shell
            L_fhead=max((2*Dnzl),(0.5*D_channel));  % entry head chamber length
        case 'C'
            D_fCoverplate=D_Flng;
            s_fChannelcover = s_channelcover;       % flat plate bolted to flanges
            L_fhead=max((2*Dnzl),(0.75*D_channel)); % entry head chamber length
        case 'N'
            D_fCoverplate=D_Flng;
            s_fChannelcover = s_channelcover;       % flat plate bolted to flanges
            L_fhead=max((2*Dnzl),(0.75*D_channel)); % entry head chamber length
    end
    switch rHead
        case '0'
            D_rCoverplate=0;                        % no cover plate
            s_rChannelcover = 0;                    % no channel
            L_rhead=0;                              % no rear head chamber 
        case 'L'
            D_rCoverplate=D_Flng;
            s_rChannelcover = s_channelcover;       % flat plate bolted to flanges
            L_rhead=0.75*D_channel;                 % entry head chamber length
        case 'M'  
            D_rCoverplate=Du_channel;
            s_rChannelcover = s_channelshell;       % bonnet cover, estimate same thickness as shell
            L_rhead=0.5*D_channel;                  % entry head chamber length
        case 'N'        
            D_rCoverplate=D_Flng;
            s_rChannelcover = s_channelcover;       % flat plate bolted to flanges
            L_rhead=0.75*D_channel;                 % entry head chamber length
        case 'T'
            D_rCoverplate=Du_channel;
            s_rChannelcover = s_channelshell;       % bonnet cover, estimate same thickness as shell
            FdimsFltHead = flangedimsFltHead((Du_shell*1000),(D_OTL*1000),(s_channel*1000),(s_tbplate*1000));
            L_rhead=0.5*D_channel;                  % entry head chamber length
            % floating head length specified later in volume calculations
        case 'U'
            D_rCoverplate=Du_channel;
            s_rChannelcover = s_channelshell;       % bonnet cover, estimate same thickness as shell
            L_rhead=0.75*D_channel;                 % entry head chamber length
    end
    Dtbsheet = Du_shell;    % Fdims.D;


    % METAL VOLUMES
    %--------------
    % volume of flanges, tube sheet
    Vflngs = Nflng*Fdims.V;                                         % volume of flanges for attaching channel cover plate and tube sheet
    Vtb_plate = s_tbplate * (0.25*pi*Dtbsheet^2);                   % tube sheet volume V [m3] BEFORE drilling, metal from drilling holes INCLUDED
    Vtb_plates = N_tbplates * Vtb_plate;
    Vbfl = N_bfls*s_bfl*AbflCond;                                    % baffle/support plate V [m3] BEFORE drilling, metal from drilling holes INCLUDED
    if (AbflCond<0)
        fprintf(1,'\n\n VITTUSAATANA!!! \n\n AbflCond = %6.4f m2 \n\n\n PERRRKELE!! ',AbflCond);
        fprintf(1,'\n\n');
    end
    Vsstrips = SStripWidth*Lshell*s_bfl * 2*HXdata.StripPairs;      % 
    if (RearHeadSmallFlange>0)
        V_FltHdFlngs = FdimsFltHead.V;
        Vflngs = Vflngs + V_FltHdFlngs;
    end


    % shell metal volume: cylinder + end cover
    switch Shell
        case 'E'
            VshellCyl = Lshell * 0.25*pi*((Du_shell^2) - ((Du_shell-s_shell)^2));	% shell metal volume V [m3]    
            VshellCover = s_shell * 0.25*pi*(Du_shell^2);
            VshellCover = VshellCover*mod((TbPasses+1),2);          % no shell cover --> zero V, if odd number of tube-side passes (2nd tubesheet at rear head, no shell cover)  
            Vshell = VshellCyl+VshellCover;
        case 'X'
            LshellCone = (2^-0.5)*((Du_shell-Du_channel)/2);
            LCylinder2 = max(Fdims.L*0.001,0.100);
            % cone in shell is part where channel-diameter shell expands to full-diameter shell to provide space for steam to flow between the bundle and the shell 
            DshellConeAvg = (Du_shell+Du_channel)/2;                
            VshellCylinder1 = (Lshell-LCylinder2-LshellCone) * 0.25*pi*((Du_shell^2) - (D_ShIn^2));	% shell metal volume V [m3]
            VshellCylinder2 = LCylinder2 * 0.25*pi*((Du_channel^2) - ((Du_channel-s_shell)^2));	% shell metal volume V [m3]
            VshellCone = LshellCone * 0.25*pi*((DshellConeAvg^2) - ((DshellConeAvg-s_shell)^2));
            VshellCover = 0.25*pi*(Du_shell^2)*s_shell;                     % shell cover metal volume [m3]
            Vshell = VshellCylinder1 + VshellCylinder2 + VshellCone + VshellCover;        
        otherwise
            fubar=1;
    end

    % metal volume of tubes, including tube length covered by baffle plates or rolled into tube plate %TotTbL = TotTbL_eff + ((s_tbplate + (s_bfl*Nbfls)) * TbPasses*NtbPerPass); % di = do-(2*s_tb);
    TotTbL = TbPasses*NtbPerPass * (LStraight+(N_tbplates*s_tbplate));
    Vtb = 0.25*pi*TotTbL*((do^2)-(di^2));

    % channel total volume 
    V_fChannelCover = s_fChannelcover * 0.25*pi*(D_fCoverplate^2);
    V_fChannelCylinder = L_fhead * 0.25*pi*((Du_channel^2) - (D_channel^2));
    if (TbPasses>1)  
        V_channeldivplate = L_fhead*D_channel*s_TbPassDivPlate;       
        if (TbPasses>2) 
            V_channeldivplate = V_channeldivplate + L_fhead*(0.5*D_channel)*s_TbPassDivPlate;
        end
    else
        V_channeldivplate = 0;
    end
    V_fChannel = V_fChannelCover+V_fChannelCylinder+V_channeldivplate;


    V_rChannelCylinder = L_rhead * 0.25*pi*((Du_channel^2) - (D_channel^2));
    switch rHead
        case 'U'
            V_rChannelCover = 0;                % no channel - shell end cover considered in shell volume calculation
            V_fltHdChannel = 0;                 % no floating head inside main rear head        
            D_rHeadDivPlate = 0;                % no pass division plates - with D=0, volume becomes 0. 
        case 'T'
            V_rChannelCover = s_rChannelcover * 0.25*pi*(D_rCoverplate^2);  % outer channel cover
            D_fltHdI = FdimsFltHead.D_fltHdI*0.001;
            D_fltHdO = FdimsFltHead.D_fltHdO*0.001;
            L_fltHd = 0.5*D_fltHdO;        
            V_fltChannelCover = s_rChannelcover * 0.25*pi*(D_fltHdO^2);
            V_fltHdChannelCylinder = L_fltHd*((D_fltHdO^2)-(D_fltHdI^2)); 
            V_fltHdChannel = V_fltChannelCover+V_fltHdChannelCylinder;
            D_rHeadDivPlate = D_fltHdI;
        otherwise
            V_rChannelCover = s_rChannelcover * 0.25*pi*(D_rCoverplate^2);   
            V_fltHdChannel = 0;
            D_rHeadDivPlate = D_rCoverplate;
    end
    V_rHeadDivplate = 0;
    if (TbPasses>2) 
        V_rHeadDivplate = L_fhead*D_rHeadDivPlate*s_TbPassDivPlate;  
        if (TbPasses>4) 
            V_rHeadDivplate = V_rHeadDivplate+L_fhead*(0.5*D_rHeadDivPlate)*s_TbPassDivPlate;  
            if (TbPasses>6) 
                V_rHeadDivplate = V_rHeadDivplate+L_fhead*(0.5*D_rHeadDivPlate)*s_TbPassDivPlate;  
            end
        end
    end
    V_rChannel = V_rChannelCover+V_rChannelCylinder+V_fltHdChannel+V_rHeadDivplate;
    Vchannels = V_fChannel+V_rChannel; %Vchannel = (Lhead * 0.25*pi*((Du_channel^2) - (D_channel^2))) + V_channelcover + V_channeldivplate;


    % remaining bits
    Vtrough = 0;                                    % trough only exists in U-tube condensers with part of the tubes submerged in condensate in a subcooling section
    Vrods = Nrods * (.25*pi*(dRod^2)) * LStraight;
    Vaircoolrbfl = 0;
    Vairvent = 0;

    % total volume before drilling (V_ini), and remaining volume after drilling (V_actual)
    V_ini = Vshell + Vtb + Vbfl + Vtb_plates  + Vchannels + Vflngs + Vtrough + Vrods; % metal volume before drilling holes
    emptyholeL = (NtbPerPass*TbPasses)*(s_tbplate*N_tbplates) + (s_bfl*N_bfls*N_HolesPerBaffle);  % metal removed by drilling
    V_actual = V_ini - (emptyholeL*0.25*pi*(do^2));

    % drilling length
    L_drill_Bfl = N_bfls*s_bfl * N_HolesPerBaffle;
    L_drill_Ts = (s_tbplate*N_tbplates) * (NtbPerPass*TbPasses);

    % external dimensions and total mass
    L_STHE1 = Lshell + s_tbplate+L_fhead+s_fChannelcover;
    if (rHead=='U')
        L_STHE = L_STHE1;
    else
        L_STHE = L_STHE1 + s_tbplate+L_rhead+s_rChannelcover;
    end
    Dtot = max(D_Flng,Du_shell);
    m_tot = V_actual*rho;
end

% complete STHX dimensions
if (fubar==1)
    result.L_tot = -1;                   
    result.D_tot = -1;         
    result.m_tot = -1;    
    
    result.m_sh = 10^9;              
    result.m_tb = 10^9;               
    result.m_bfl = 10^9;               
    result.m_ts = 10^9;           
    result.m_ch = 10^9;           
    result.m_flng = 10^9;              
    result.m_other = (10^9)*rho;   
    result.s_bfl = 10^6;
    result.s_tbplate = 10^6;
    result.s_fChannelcover = 10^6;
    result.s_rChannelcover = 10^6;
    result.s_channelshell = 10^6;
    result.s_shell = 10^6;
    result.L_shell = 10^6;
    result.L_tbTot = 10^12;
    result.L_fhead = 10^6;
    result.L_rhead = 10^6;
    result.s_TbPassDivPlate = 10^6;
    result.L_drill_Bfl = 10^9;
    result.L_drill_Ts = 10^9;     
    %result.L_cut_acBfl = 2*(0.5*pi*D_aircoolrbfl+(1.1*s_bfl));        
    result.D_ShO = 10^6;            
    result.D_ShIn = 10^6;             
    result.D_Flng = 10^6;      
    result.D_ChIn = 10^6;  
    result.D_ChO = 10^6;
else
    result.L_tot = L_STHE;                   
    result.D_tot = Dtot;         
    result.m_tot = m_tot;

    % masses of individual components
    result.m_sh = Vshell*rho;               % shell cylinder + cover
    result.m_tb = Vtb*rho;                  % mass of heat transfer tubes
    result.m_bfl = Vbfl*rho;                % mass of baffle/support plates (before drilling)
    result.m_ts = Vtb_plates*rho;           % mass tube sheet (before drilling)
    result.m_ch = Vchannels*rho;            % mass of channel: cylinder and channel cover
    result.m_flng = Vflngs*rho;             % mass of flanges; 2 for channel, 1 for shell in U-tube
    result.m_other = (Vtrough+Vrods+Vairvent+Vaircoolrbfl+Vsstrips)*rho;   % spacer/tie rods, air venting tube, a possible trough for subcooling

    % dimensions of selected components
    result.s_bfl = s_bfl;
    result.s_tbplate = s_tbplate;
    result.s_fChannelcover = s_fChannelcover;
    result.s_rChannelcover = s_rChannelcover;
    result.s_channelshell = s_channelshell;
    result.s_shell = s_shell;
    result.L_shell = Lshell;
    result.L_tbTot = TotTbL;
    result.L_fhead = L_fhead;
    result.L_rhead = L_rhead;
    result.s_TbPassDivPlate = s_TbPassDivPlate;
    result.L_drill_Bfl = L_drill_Bfl;
    result.L_drill_Ts = L_drill_Ts;     
    %result.L_cut_acBfl = 2*(0.5*pi*D_aircoolrbfl+(1.1*s_bfl));        
    result.D_ShO = Du_shell;            
    result.D_ShIn = D_ShIn;             
    result.D_Flng = D_Flng;      
    result.D_ChIn = D_channel;  
    result.D_ChO = Du_channel;
end

%tubes
result.LStraight = HXdata.LStraight;
result.do = do;                         % [m]
    
result.fock = fubar;
result.StripPairs = HXdata.StripPairs;
result.SStripWidth = HXdata.SStripWidth;
result.NtbPerPass = HXdata.NtbPerPass;


if(plotti>0)
    fid = fopen(outputfname,'w');     
    fseek(fid, 0, 'eof');
    fprintf(fid,'\n\n');
    fprintf(fid,'tube outside diametre [mm]:         %4.1f\n',do*1000);
    fprintf(fid,'tube wall thickness [mm]:           %3.1f\n',s_tb*1000);
    fprintf(fid,'tube pitch [mm]:                    %4.1f\n\n',P);

    fprintf(fid,'Baffle plates:                      %1.0f\n',N_bfls);
    %fprintf(fid,'Baffle plates, subcooling section: %1.0f\n\n',Nbfl_A);

    fprintf(fid,'baffle plate thickness [mm]:        %4.1f\n',s_bfl*1000);
    fprintf(fid,'tube plate thickness [mm]:          %5.1f\n',s_tbplate*1000);
    fprintf(fid,'shell thickness [mm]:               %4.1f\n',s_shell*1000);
    fprintf(fid,'front head end plate thickness [mm]:%4.1f\n',s_fChannelcover*1000);
    fprintf(fid,'rear head end plate thickness [mm]: %4.1f\n',s_rChannelcover*1000);
    fprintf(fid,'channel shell thickness [mm]:       %5.1f\n',s_channelshell*1000);
    fprintf(fid,'front head channel length [mm]:     %4.0f\n\n',L_fhead*1000);
    fprintf(fid,'rear head channel length [mm]:      %4.0f\n\n',L_rhead*1000);

    fprintf(fid,'shell length [mm]:                %5.0f\n',Lshell*1000);
    fprintf(fid,'total length [mm]:                %5.0f\n',L_STHE*1000);
    fprintf(fid,'total diametre [mm]:              %4.0f\n',Dtot*1000);
    fprintf(fid,'total mass [kg]:                  %6.0f\n',m_tot);
    fclose(fid);
end

end


function s_plate = platethickness(p_tube,p_shell,Ds,Dci,s_sh,s_ch,do,P,metalprops)

% p_tube = max tube-side pressure (+ or -) [MPa]
% p_shell = max shell-side pressure (+ or -) [MPa]
% ds = shell diameter, inside [mm]
% l = shell length [mm]
% do = tube outside diameter [mm]
% P = tube pitch [mm]
% a1 =  shell thickness lower bound [mm]
% b1 =  shell thickness upper bound [mm]

sigmal=metalprops.sigmal;
n=metalprops.n;

c1 = 2.0;                                   % korroosiovara
c2 = 0.5;                                   % valmistustoleranssi
%n=1.5;                                     % varmuuskerroin

% TUBE PLATE THICKNESS
%---------------------
vR = (P-do)/P;                              % reikäkentän vaikutus
D = Ds + s_sh;                              % laskentahalkaisija
%C = 0.40;                                  % tube plate between flanges; SFS 2774 Taulukko 2 (s.15)
C = 0.50;                                   % tube plate welded from one side to the flange, s_ts>3s_sh; SFS 2774 Taulukko 2 (s.14)
vh = 1;                                     % welds non-destructively inspected -> wel strength factor can be assumed 1 
v = vh*vR;

p=max(p_tube,p_shell);
s_tbplate = C*D * ((p*n / (v*sigmal) )^0.5);
s_tbplate = s_tbplate + c1 + c2;


% CHANNEL COVER PLATE THICKNESS
%------------------------------
C=.35;                                      % flange attachment, SFS 2615 8.4 Kuva 15 (s.11)
v=1;                                        % no welds, no holes

Dco = Dci + s_ch; 
s_cvrplate = C*Dco * ((p*n / (v*sigmal) )^0.5);
s_cvrplate = s_cvrplate + c1 + c2;

s_plate(1)=s_tbplate;
s_plate(2)=s_cvrplate;
end


function fx = flangedimsSimple(D,s,s_ts)
% D = shell inside diameter [mm]
% s = shell wall thickness [mm]
% s_ts = tubesheet thickness [mm]
% 

Du = D+(2*s);
D4 = Du + 2.5*s_ts;
D3 = Du + 0.5*s_ts;
L_hub = 2*s_ts;
%D3 = (1.04*D) + 21;
%D4 = (1.10*D) + 76;
%L = min(((76*log(D)) - 399), (s_ts));
t = (56*log(D)) - 304;
t = min(t,s_ts);
L = t+L_hub;

%L_hub = L-t;
D_hub_avg = 0.5*(D3+Du);

Vflng_1 = t * ((D4^2)-(Du^2));
Vflng_hub = L_hub * ((D_hub_avg^2)-(Du^2));

V = Vflng_1 + Vflng_hub;                    % total volume of single flange [mm3]
V = V*(.001^3);                             % total volume of single flange [m3]

fx.V=V;
fx.D=D4*0.001;                              % channel cover plate outside diameter [mm]
fx.L=L;                                     % flange "foot" length
end

function fx = flangedimsFltHead(D_sh,D_OTL,s,s_ts)
% D_sh = shell inside diameter [mm]
% D_OTL = shell inside diameter [mm]
% s = channel wall thickness [mm]
% s_ts = tubesheet thickness [mm]

Du = D_OTL+0.020+(2*s);
D4 = D_sh;
D3 = 0.65*Du + 0.35*D4;
L_hub = s_ts;
t = (56*log(Du)) - 304;
t = min(t,s_ts);
L = t+L_hub;

%L_hub = L-t;
D_hub_avg = 0.5*(D3+Du);

Vflng_1 = t * ((D4^2)-(Du^2));
Vflng_hub = L_hub * ((D_hub_avg^2)-(Du^2));

V = Vflng_1 + Vflng_hub;                    % total volume of single flange [mm3]
V = V*(.001^3);                             % total volume of single flange [m3]

fx.V=V;
fx.D=D4*0.001;                              % channel cover plate outside diameter [mm]
fx.L=L;                                     % flange "foot" length
fx.D_fltHdO = Du;
fx.D_fltHdI = Du-(2*s);
end


function s_sh = channelshell_intp(p,ds,metalprops)
% calls s_internalp(): SFS 2611, 1996:
% p = overpressure inside [MPa]
% ds = shell diameter [mm]
% sigma = laskentalujuus [N/mm2]
% s_sh = shell minimum wall thickness [mm]

sigmal = metalprops.sigmal;
n = metalprops.n;

if (isfield(metalprops,'c1'))
    c1=metalprops.c1;                       % [mm] corrosion allowance
else
    c1 = 3.0;                            	% [mm] korroosiovara
end

if (isfield(metalprops,'c2'))
    c2=metalprops.c2;                       % [mm] manufacturing tolerance
else
    c2 = 2.0;                               % [mm] työstötoleranssi
end

s_sh = s_internalp(p,ds,sigmal,n) + c1 + c2;

end

function s_sh = shellthickness(p_out,p_in,ds,l,s_a1,s_b1,metalprops)

% p_out = max overpressure outside [MPa]
% p_in = max overpressure inside [MPa]
% ds = shell diameter, inside [mm]
% l = shell length [mm]
% s_a1 =  shell thickness lower bound [mm]
% s_b1 =  shell thickness upper bound [mm]

if (isfield(metalprops,'c1'))
    c1 = metalprops.c1;                     % corrosion allowance
else
    c1 = 3.0;                            	% korroosiovara
end

if (isfield(metalprops,'c2'))
    c2 = metalprops.c2;                     % manufacturing tolerance
else
    c2 = 2.0;                               % työstötoleranssi
end

Et = metalprops.Et;                         % N/mm2
nyy = metalprops.nyy;                       % [-]
nk = metalprops.nk;                         % [-] nk=3; SFS 2862, 5.1 Kimmoinen lommahdus (pp.184)
sigmal = metalprops.sigmal;                 % N/mm2
n =  metalprops.n;                          % [-] 

if (p_out==0)
    s_sh_a(1)=0;
    s_sh_a(2)=0;
else    
    s_sh_a(1) = s_elasticbuckling(p_out,ds,l, s_a1, s_b1, Et,nyy,nk);
    s_sh_a(2) = s_nonelastickbuckling(p_out,ds,l, s_a1, s_b1, sigmal,n);
end
s_sh_a(3) = s_internalp(p_in,ds,sigmal,n);  % SFS 2611 

s_sh = max(s_sh_a) + c1 + c2;       

end





function s_sh = s_internalp(p,ds,sigmal,n)
% SFS 2611, 1996:
% p = overpressure inside [MPa]
% ds = shell diameter [mm]
% sigma = laskentalujuus [N/mm2]
% s_sh = shell minimum wall thickness [mm]

s_sh = ds*p / ( (2*sigmal/n) - p);

end


function s_sh = s_nonelastickbuckling(p,ds,l, s_a1, s_b1, sigmal,n)
% out-of-roundness assumed at e = 0.015;
% jos suurempi epäpyöreys, ks. SFS 2862 pp.186, alalaita --> kuva 7.

if (ds>5*l)                                     % SFS 2862 kpl 8.2, erikoistapaus hyvin lyhyille / for very short  only
    s_sh_a(1) = ds / (2*((sigmal/p)-1));        % SFS 2862 kpl 8.2 kimmoton lommahdus, yht.(7)
    s_sh_a(2) = l * (((p*n)/(3*sigmal))^0.5);   % SFS 2862 kpl 8.2 kimmoton lommahdus, yht.(8)
    s_sh = min(s_sh_a);                         % SFS 2862 s.4 (186), kpl 8.2.1: "ei saa alittaa PIENEMPÄÄ kaavoista (7),(8)

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  GOLDEN SECTION SEARCH                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = (3 - (5^0.5)) / 2;

    s_a = s_a1;
    s_b = s_b1;
    s_c = s_a + gamma*(s_b - s_a);
    s_d = s_b - gamma*(s_b - s_a);

    p_c = p_nonelasticbuckling(n,sigmal,s_c,(ds+2*s_c),l);
    f_c = abs(p-p_c);

    p_d = p_nonelasticbuckling(n,sigmal,s_d,(ds+2*s_c),l);
    f_d = abs(p-p_d);

    continuing = 1;

    fevals=2;
    while continuing == 1
        if (f_c < f_d)
            s_b = s_d;

            s_d = s_c;
            p_d = p_c;
            f_d = f_c;        

            s_c = s_a + gamma*(s_b-s_a);
            p_c = p_nonelasticbuckling(n,sigmal,s_c,(ds+2*s_c),l);
            f_c = abs(p-p_c);
        else
            s_a = s_c;

            s_c = s_d;
            p_c = p_d;
            f_c = f_d;

            s_d = s_b - gamma*(s_b-s_a);
            p_d = p_nonelasticbuckling(n,sigmal,s_d,(ds+2*s_c),l);
            f_d = abs(p-p_d);
        end
        fevals = fevals+1;
        if ((.01*round(100*s_a)) == (.01*round(100*s_b)))
            continuing=0;
        end
    end

    s_sh=0.5*(s_a+s_b);    
end
end

function p = p_nonelasticbuckling(n,sigmal,s,du,l)

% SFS 2862: Paineastiat, mitoitus. SFS-käsikirja 13, pp.186, eq.(9)

sperdu = s/du;
aux = (1 - (du/(5*l))) / sperdu;
p = (2*sigmal/n) * sperdu * 1/(1 + 1.5*0.015*aux);

end


function s_sh = s_elasticbuckling(p,ds,l, s_a1, s_b1, Et,nyy,nk)

% p = max overpressure outside [MPa]
% ds = shell diameter, inside [mm]
% l = shell length [mm]
% s_a1 =  shell thickness lower bound [mm]
% s_b1 =  shell thickness upper bound [mm]
% Et = Young's modulus [N/mm2]
% nyy = Poisson's ratio [-]
% nk = safety factor [-]

%s = 15;
C=pi*ds/(2*l);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  GOLDEN SECTION SEARCH                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = (3 - (5^0.5)) / 2;

s_a = s_a1;
s_b = s_b1;
s_c = s_a + gamma*(s_b - s_a);
s_d = s_b - gamma*(s_b - s_a);

sol_c = Nmin(Et,nk,C,s_c,ds,l,nyy);
p_c = sol_c(1);
f_c = abs(p-p_c);

sol_d = Nmin(Et,nk,C,s_d,ds,l,nyy);
p_d = sol_d(1);
f_d = abs(p-p_d);

continuing = 1;

fevals=2;
while continuing == 1
    if (f_c < f_d)
        s_b = s_d;
        
        s_d = s_c;
        p_d = p_c;
        f_d = f_c;        
        
        s_c = s_a + gamma*(s_b-s_a);
        sol_c = Nmin(Et,nk,C,s_c,ds,l,nyy);
        p_c = sol_c(1);
        f_c = abs(p-p_c);
    else
        s_a = s_c;
        
        s_c = s_d;
        p_c = p_d;
        f_c = f_d;
        
        s_d = s_b - gamma*(s_b-s_a);
        sol_d = Nmin(Et,nk,C,s_d,ds,l,nyy);
        p_d = sol_d(1);
        f_d = abs(p-p_d);
    end
    fevals = fevals+1;
    if ((.01*round(100*s_a)) == (.01*round(100*s_b)))
        continuing=0;
    end
end

s_sh=0.5*(s_a+s_b);
end

function sol = Nmin(Et,nk,C,s_trial,ds,l,nyy)
    % enumerate through N=1,...,10 buckling waves
    %N = [2 3 4 5 6 7 8 9 10];
    %p_N = zeros(1,9);
    N0 = 2+ceil(1.63*sqrt((ds/l)*sqrt(ds/s_trial)));
    N = linspace(1,N0,N0);
%     for i = 1:4
%        N(i)=floor(N0)-2+i;
%     end
    p_N = zeros(1,N0);
    Nmin = 99;
    pmin = 9999;
    for i = 1:N0
        p_N(i) = p_elasticbuckling(Et,nk,N(i),C,s_trial,((2*s_trial)+ds),nyy); % [MPa]
        if p_N(i) < pmin
            pmin = p_N(i);
            Nmin = i;
        end
    end
    %plot(N(2:N0),p_N(2:N0))
    sol(1) = pmin;
    sol(2) = Nmin;
end

function p = p_elasticbuckling(Et,nk,N,C,s,du,nyy)

N21 = (N^2)-1;
sperdu = s/du;
f_ny = 2 / (3*(1-(nyy^2)));

a = 2 / (N21*((1 + ((N/C)^2))^2));
b = N21 + (((2*(N^2)) - 1 - nyy) / (1 + ((N/C)^2) ));

p = (Et/nk) * ((a*sperdu) + (f_ny*b*(sperdu^3)));
end