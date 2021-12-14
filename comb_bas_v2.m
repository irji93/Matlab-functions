  function [MFB_m,Qb_m,p_m,RGF_out,ut_m,Sl_m]=comb_bas_v2(C,in)    
% Load info
ang=-180:0.2:540-0.2;                                                       % Angle discretization
D=72e-3;                                                                    % Bore (m)
Lb=128e-3;                                                                  % Connecting-rod length (m)
Lm=40.6e-3;                                                                 % Crankshaft radius (m)
Kc=10.6;                                                                    % Compression ratio (-)

Vdis=pi*D^2/4*Lm*2;                                                         % Displaced volume (m^3)
Vcc=Vdis/(Kc-1);                                                            % Minimum volume (dead volume) (m^3)
x1=Lb+Lm-Lm*cos(ang*pi/180)-sqrt(Lb^2-Lm^2*sin(ang*pi/180).^2);             % Piston position throughout the cycle (angle in rad)
hc=(x1+Vcc/(pi*D^2/4))';                                                    % Chamber height throughout the cycle (m)
V=pi*D^2/4*x1'+Vcc;                                                         % Chamber volume throughout the cycle (m^3)
dV=filter([3 2 1 0 -1 -2 -3],28,[V(end-2:end);V;V(1:3)]);dV=dV(7:end);      % Use the rational transfer function to create dV vector

me_m=zeros(size(V)); 
mb_m=zeros(size(V));
Qb_m=zeros(size(V)); 
dQw_m=zeros(size(V)); 
rb_m=zeros(size(V)); 
Af_m=zeros(size(V));
Vb_m=zeros(size(V));
T=zeros(size(V));

% Rename inputs
mair=in.mair;
mfuel=in.mfuel;
SAn=in.SAn;
p0=in.p0;
lambda=in.lambda;
n=in.n;
RGF_m=in.RGF;
% IVC=in.IVC;                                                                 % Intake Valve Closing (º)
EVOn=in.EVOn;                                                                 % Exhaust Valve Opening (º)
EVCn=in.EVCn;                                                                 % Exhaust Valve Closing (º)
IVCn=in.IVCn;                                               % Intake Valve Closing Position
                                              % Exhaust Valve Closing Position

% C=in.C;
ce_m=in.ce;                                                                           

% Load variability
C1=C(1);                                                              % Sl variability
C2=C(2);                                                             % ut variability 
C3=C(3);                                                           % lm variability
rb_m(SAn)=5e-4;%5e-4                                                            % Initial Radio variability

% Model 
it=SAn;                                                                     % Combustion start position
Fs=6*5*n;                                                                   % Sampling frecuency (at constant rpm)
dt=1./Fs;                                                                   % Time step (Time between samples)
kappa=in.kcomp;                                                                  % Polytropic coefficient
pm_m=p0.*(V(IVCn)./V).^kappa;
p_m=p0.*(V(IVCn)./V).^kappa;

     
m_m=(mair+mfuel);
Fst=mair./mfuel;% Mixture mass (mfuel + mair).
mub_m=ones(size(V))*m_m;                                                    % Unburned mass
Vub_m=V;                                                                    % Unburned volume
dub_m=m_m./Vub_m;                                                           % Unburned density
Tub_m=p_m./dub_m/287;                                                       % Unburned temperature
vis_m=3.3e-7*Tub_m.^0.7;  % Dynamic viscosity 

% whoscni parameters
Cw1=2.28;Cw2=3.24e-3;
K =Vdis*Tub_m(IVCn)./(p_m(IVCn)*V(IVCn));
Sp = 2*Lm*n/60;Ap = pi*D^2/4;Tw=100;

%model paremetrs
a1=(2.4-0.271*lambda^3.51);                                                 % Alpha
a2=(-0.357+0.14*lambda^2.77);                                               % Beta
a3=C1*(30.5-54.9.*(lambda-1.21)^2)/100;                                     % Reference laminar flame speed (Sl0) (m/s)
Sl_m=a3.*(Tub_m./298).^a1.*(p_m(:)./101e3).^a2.*(1-2.06.*RGF_m.^0.77);      % Laminar flame speed power law model with residual gas effects (m/s)
up=4*Lm.*n/60;                                                              % Mean piston velocity and adjustable constant
ut0=C2*up*(dub_m(SAn)/dub_m(IVCn))^(1/2);                                   % Turbulent intensity at start of combustion
ut_m=ut0*(dub_m./dub_m(SAn)).^(1/3);                                        % Turbulent intensity
L0=C3^2*hc(SAn);                                                            % Initial length scale
L_m=L0.*(dub_m(SAn)./dub_m).^(1/3);                                         % Characteristic length scale
lm_m=sqrt(vis_m.*L_m./ut_m./dub_m);                                         % Taylor microscale 

%% First instant after the spark with the selected burning radius        
Vb_m(SAn)=2/3*pi*rb_m(SAn)^3;                                               % Recalculated burned volume for next time step
Vub_m(SAn)=V(SAn)-Vb_m(SAn);                                                % Recalculated unburned volume for next time step
mb_m(SAn)=m_m*(1-Vub_m(SAn)/V(SAn))/(3*Vub_m(SAn)/V(SAn)+1);                % Recalculated burned mass for next time step
mub_m(SAn)=m_m-mb_m(SAn);                                                   % Recalculated unburned mass for next time step  
dub_m(SAn)=mub_m(SAn)/Vub_m(SAn);                                           % Recalculated burned gas density for next time step         
% Qb_m(SAn)=(mb_m(SAn)-mb_m(SAn-1))/(1+14.6*lambda)*(1-RGF_m)*46e6*ce_m;      % Recalculated energy release in each time step (J/samp)
% p_m(SAn)=p_m(SAn-1)+(Qb_m(SAn)-p_m(SAn-1).*dV(SAn-1)*kappa/(kappa-1))./V(SAn-1)*(kappa-1); % Recalculated pressure for next time step
Tub_m(SAn)=p_m(SAn)./dub_m(SAn)/287;                                        % Recalculated unburned gas temperature for next time step
Af_m(SAn)=2*pi*(rb_m(SAn))^2;                                               % Superficial flame area for next time step
vis_m(SAn)=3.3e-7*Tub_m(SAn).^0.7;                                          % Dynamic viscosity for next time step
Sl_m(SAn)=a3.*(Tub_m(SAn)./298).^a1.*(p_m(SAn)./101e3).^a2.*(1-2.06.*RGF_m.^0.77); % Laminar flame speed for next time step (m/s) 
ut_m(SAn)=ut0*(dub_m(SAn)./dub_m(SAn)).^(1/3);                              % Turbulence intensity for next time step
L_m(SAn)=L0.*(dub_m(SAn)./dub_m(SAn)).^(1/3);                               % Characteristic length scale for next time step
lm_m(SAn)=sqrt(vis_m(SAn).*L_m(SAn)./ut_m(SAn)./dub_m(SAn));                % Taylor microscale for next time step

%% Loop while the unburned mass gets burned 
while it<=EVOn+50   
    % Computation of the entrained mass
    if me_m(it)+dt*(dub_m(it).*Af_m(it).*(ut_m(it)+Sl_m(it)))<m_m       % If entrained mass < total mixture mass
        me_m(it+1)=me_m(it)+dt*(dub_m(it).*Af_m(it).*(ut_m(it)+Sl_m(it))); % Entrained mass in the next time step     
    else
        me_m(it+1)=m_m;                                                 % All the mixture mass has been entrained
    end

    % Computation of the burned mass
    if mb_m(it)+dt*((me_m(it)-mb_m(it))/lm_m(it)*Sl_m(it)+dub_m(it)*Af_m(it)*Sl_m(it))<=m_m % If there is still entrained mass not burned
        mb_m(it+1)=mb_m(it)+dt*((me_m(it)-mb_m(it))/lm_m(it)*Sl_m(it)+dub_m(it)*Af_m(it)*Sl_m(it));% Burned mass in the next time step
    else
        mb_m(it+1)=m_m;                                                 % If the burned mass is bigger than the entrained mass
    end

    % Next iteration update variables
    Qb_m(it+1)=(mb_m(it+1)-mb_m(it))/(Fst+1)*45e6*ce_m; % Energy release in each time step (J/samp) with combustion efficiency (0.7)
    T(it+1)=p_m(it).*V(it)./286/m_m;
    g=1.38-0.2*exp(-900./T(it));
    % whoscni parametrs
    w= Cw1*Sp+Cw2*K*(p_m(it+1)-pm_m(it+1));    
    h1 = 0.013*D^-0.2*p_m(it+1).^0.8.*T(it+1).^-0.53.*w.^0.8;
    A = 2*Ap+pi*D*hc(it+1);  
    dQw_m(it+1)= A/(6*n).*h1.*(T(it+1)-Tw)/5*0.75;        
    p_m(it+1)=p_m(it)+((Qb_m(it+1)-dQw_m(it+1))-p_m(it).*dV(it)*g/(g-1))./V(it)*(g-1); % Next pressure value from heat equation
    mub_m(it+1)=m_m-mb_m(it+1);                                         % Unburned mass for next time step
    dub_m(it+1)=(m_m+3*mb_m(it+1))./V(it+1);                            % According to Heywood
    Vub_m(it+1)=mub_m(it+1)/dub_m(it+1);                                % Unburned volume for next time step                              
    Tub_m(it+1)=p_m(it+1)./dub_m(it+1)/287;                             % Next unburned gas temperature from ideal gas equation 
    vis_m(it+1)=3.3e-7*Tub_m(it+1).^0.7;                                % Dynamic viscosity for next time step
    Vb_m(it+1)=V(it+1)-Vub_m(it+1);                                     % Burned gas volume for next time step  
    rb_m(it+1)=real((Vb_m(it+1)*3/2/pi).^(1/3));                        % Burned gas radius for next time step
    if rb_m(it+1) > hc(it+1)                                            % Loop for semi-sphere radius > chamber height
        rb_m(it+1)=sqrt((Vb_m(it+1)+pi*hc(it+1)^3/3)/pi/hc(it+1));      % Recalculated burning radius for next time step ; 
        Af_m(it+1)=2*pi*rb_m(it+1)*hc(it+1);                            % Superficial flame area for next time step
    else
        Af_m(it+1)=2*pi*rb_m(it+1)^2;                                   % Superficial flame area for next time step
    end

    % Laminar flame speed
    Sl_m(it+1)=a3.*(Tub_m(it+1)./298).^a1.*(p_m(it+1)./101e3).^a2.*(1-2.06.*RGF_m.^0.77); % Laminar flame speed for next time step (m/s) 
    % Turbulent intensity
    ut_m(it+1)=ut0*(dub_m(it+1)./dub_m(SAn)).^(1/3);                    % Turbulence intensity for next time step
    % Microscale length
    L_m(it+1)=L0.*(dub_m(SAn)./dub_m(it+1)).^(1/3);                     % Characteristic length scale for next time step
    lm_m(it+1)=sqrt(vis_m(it+1).*L_m(it+1)./ut_m(it+1)./dub_m(it+1));   % Taylor microscale for next time step

    it=it+1;                                                            % Next step      
end
EVOni = EVOn+50;
SPO = EVOni+300;

p_m(EVOni:SPO)=p_m(EVOni)+(p_m(SPO)-p_m(EVOni))./(ang(SPO)-ang(EVOni)).*(ang(EVOni:SPO)-ang(EVOni));
p_m(SPO:end) = p0;



RGF_out=(p_m(EVCn)./p_m(EVOn)).^(1/kappa)*(V(EVCn)/V(EVOn));
MFB_m=cumsum(Qb_m(in.c1:in.c2))/sum(Qb_m(in.c1:in.c2));
end