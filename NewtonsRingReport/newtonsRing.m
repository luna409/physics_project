close all; 
figure('Position',[90 164 873 483]); 

%change the value of L,R,H
L=580; %wave length
R=5; %curvature radius
H=0; %Air thickness

a1=axes('Position',[0.4,0.16,0.4,0.7]); 
[x,y]=meshgrid(linspace(-0.005,0.005,200)); 
r2=(x.^2+y.^2); 
Di=(2*H+2*(R-sqrt(R^2-r2))*1e9)/L;
In=abs(cos(Di*pi*2));
if (380<=L) && (L<440)
    cr = -1.*(L-440.)/(440.-380.);
    cg = 0.;
    cb = 1.;
elseif (440<=L) && (L<490)
        cr = 0.;
        cg= (L-440.)/(490.-440.);
        cb= 1;
elseif (490<=L) && (L<510)
        cr = 0.;
        cg= 1.;
        cb= -1.*(L-510.)/(510.-490.);
elseif (510<=L) && (L<580)
        cr = (L-510.)/(580.-510.);
        cg= 1.;
        cb= 0.;
elseif (580<=L) && (L<645)
        cr = 1.;
        cg = -1.*(L-645.)/(645.-580.);
        cb = 0.;
else
    (645<=L) && (L<780)
        cr = 1.;
        cg = 0.;
        cb = 0.;
end

Ik(:,:,1)=In*cr;
Ik(:,:,2)=In*cg; 
Ik(:,:,3)=In*cb; 
Pc=imshow(Ik,[]); 
title('Simulate Ring','fontsize',18); 
 
TT=In*cr+In*cg+In*cb;

%wavelength title
Lt=uicontrol(gcf,'style','text',... 
'unit','normalized','position',[0.06,0.8,0.23,0.06],... 
'BackgroundColor',0.7*[1,1,1],'ForegroundColor',[1,1,1],... 
'string',sprintf('Wavelengh:%3.1f nm',L),'fontsize',16,'fontname','timesnewroman');


%Curvature radius title
Rt=uicontrol(gcf,'style','text',... 
'unit','normalized','position',[0.06,0.6,0.23,0.06],... 
'BackgroundColor',0.7*[1,1,1],'ForegroundColor',[1,1,1],... 
'string',sprintf('Curvature radius:%1.1f m',R),'fontsize',16,'fontname','timesnewroman'); 


%Air thickness title
Ht=uicontrol(gcf,'style','text',... 
'unit','normalized','position',[0.06,0.4,0.23,0.06],... 
'BackgroundColor',0.7*[1,1,1],'ForegroundColor',[1,1,1],... 
'string',sprintf('Air thickness:%1.1f nm',H),'fontsize',16,'fontname','timesnewroman'); 

