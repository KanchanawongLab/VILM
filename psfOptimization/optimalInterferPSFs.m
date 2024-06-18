% change the folder to current one to run the script

clear all
close all

global addvortex %add vortex means whether to optimize with vortex phase mask
addvortex=0;     % default 0, geenerate optz psf, 1 is set to generate vxz psf.
global showPSFs; % whether to show PSFs
showPSFs=0;      % default 0, not display, 1 display psfs and crlb.

global kxx1 kyy1 kz_c1 kz_c2
global CTF1
global R L NN MM 
global zpos
global zern_p
global AA bg
global k0;


%% basic parameters of the imaging system
na=1.33;  mag=60;lambda=0.670; 
%z=(0:0.1:1); % z's dimension should be 1*num_z1 
d=50; d1=[]; d2=[];
pix=7.2;  spix=pix/mag;
fpupil=[]; nn1=5; 
spix=spix/1;  nn1=nn1*1;
num_x=2*nn1+1;
x0=0;y0=0; n1=1.52; n0=1.33;
k0=2*pi/lambda; kmax=k0*na;


%% PSF is calculated using scalarczt method, below are setting parameters for it
[kx_c xx_c]=czt_getcztcoord(num_x,300);
kx=(k0*na).*kx_c; 
xx=(num_x*spix/2).*xx_c;
[kxx1 kyy1]=meshgrid(kx);
[pxx1 pyy1]=meshgrid(xx);
NN=length(kx); MM=length(xx);
R=(k0*na) ;   L=(num_x*spix/2);
CTF1=circ(kxx1,kyy1,2*k0*na); 
%kz in the immersion oil
kz_c1=CTF1.*sqrt((k0*n1)^2-kxx1.^2-kyy1.^2);
%kz in the medium
kz_c2=CTF1.*sqrt((k0*n0)^2-kxx1.^2-kyy1.^2);
do_apod=1;
if do_apod==1
apod=sqrt(kz_c1/(n1.*k0));
else 
    apod=ones(nn,nn);
end

%% range -300nm to 300nm
dz=0.01;
zpos=(-30:30).*dz;


%% introduce the phase zernike
R=(k0*na) ; 
pupil_basic=CTF1;
% generate the zernike base
[theta,r] = cart2pol(kxx1./R,kyy1./R);
idx = r<=1;
% zern_p=[];
% p=3:35 as the 1,2 means translate
for p=3:35
    temp_zern1= zeros(size(kxx1));
    % z = nan(size(kxx));
    temp_zern = zernfun2(p,r(idx),theta(idx));
    %disp(sum(temp_zern(:)));
    temp_zern1(idx)=temp_zern;
    zern_p(:,:,p-2) = temp_zern1;
    %showing the zernike pattern
    figure(11)
    imagesc(CTF1.*zern_p(:,:,p-2));
    colormap('hot');
end
      
%setting the optimization 
zern_no_short=28-3+1;
%choose random or zeros as the initial guess.
fit_par_short(1:zern_no_short)=0.*rand(1,zern_no_short);
obj_ipalmpsfoptim(0.*fit_par_short);

lb(1:zern_no_short)=-3.*ones(1,zern_no_short); 
ub(1:zern_no_short)= 3.*ones(1,zern_no_short);
   
%initialization using simulated annealing     
options_psfoptim = optimoptions(@fmincon,'MaxIterations',100);
Afmin=[];bfmin=[];
Aeqfmin=[];beqfmin=[];
fit_par_short=fmincon(@obj_ipalmpsfoptim,fit_par_short,Afmin,bfmin,Aeqfmin,beqfmin,lb,ub,[],options_psfoptim);   











