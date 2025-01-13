% change the folder to current one to run the script

clear all;  
close all;

% parameters 
% na: numerical aperture,  lambd: wavelength of the laser
% spixel: pixel size for excitation psf  
% num_x : size of the excitation psf
% LL 
% PSFsmooth:whether to generate pixel smoothed PSFs, default no smooth with
% value 0

%% basic parameters for calculating the excitation point spread function
na=1.33;      % numerical aperture of the optical system  
lambda=0.67;  % micrometer
%z=(0:0.1:1); % z's dimension should be 1*num_z1 
d=50; d1=[]; d2=[];
pix=7.2;%  micrometer  
mag=60;
spix=pix/mag;  % pixel size 
num_x=2*10+1; % number of pixels in the excitation point spread function 
PSFsmooth=0;  % change it to 1 to generate smoothed PSFs

%whether generate smoothed PSfs
if PSFsmooth==1
    spix=spix/10;
    num_x=2*100+1;
end

%% define the coordinats in the fourier plane and sample plane
n1=1.52; n0=1.33;   %  refractive index for oil and sample
k0=2*pi/lambda; fpupil=[];
kmax=k0*na;
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
% define the apodization function for calculating the excitation point
% spread function

do_apod=1;
if do_apod==1
kz_c11=sqrt((k0*n1)^2-kxx1.^2-kyy1.^2);
apod=1./sqrt(kz_c11/(n1.*k0));
else 
    apod=ones(NN,NN);
end

%% calculating the PSFsz range -150nm to 150nm
dz=0.003;
zpos=(-50:50).*dz;

ddxy=1e-4; % 0.1nm
ddz=1e-4;  %0.1nm
x0=0; y0=0;
z=0;
z1=0;z2=0;
fai=atan2(kxx1,kyy1);
for jj=1:length(zpos)
    focpos=zpos(jj);
    %focpos=0;
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z,focpos);
    fmask1=exp(1i.*fai);  
    fkxky=fpupil.*fmask1.*apod;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z1,z2-focpos);
    fmask1=exp(-1i.*fai);  
    fkxky=fpupil.*fmask1.*apod;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsf(:,:,jj,1)=abs(im+im1).^2;
    beadpsf(:,:,jj,2)=abs(im-im1).^2;

    h1=figure(1);
    set(h1,'Position',[100 100 600 300])
    subplot(121)
    imagesc(abs(im+im1).^2);
    title('PSF in channel one')
    colormap('copper')
    subplot(122)
    imagesc(abs(im-im1).^2);
    title('PSF in channel two')
    colormap('copper')
    


    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+ddxy,y0,z,focpos);
    fmask1=exp(1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+ddxy,y0,z1,z2-focpos);
    fmask1=exp(-1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsfdx(:,:,jj,1)=abs(im+im1).^2;
    beadpsfdx(:,:,jj,2)=abs(im-im1).^2;


    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0+ddxy,z,focpos);
    fmask1=exp(1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0+ddxy,z1,z2-focpos);
    fmask1=exp(-1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsfdy(:,:,jj,1)=abs(im+im1).^2;
    beadpsfdy(:,:,jj,2)=abs(im-im1).^2;


    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z,focpos+ddz);
    fmask1=exp(1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z1,z2-focpos-ddz);
    fmask1=exp(-1i.*fai); 
    fkxky=fpupil.*fmask1.*apod;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsfdz(:,:,jj,1)=abs(im+im1).^2;
    beadpsfdz(:,:,jj,2)=abs(im-im1).^2;


end


%%normalizd the PSF, multiply the intensity add, the noise
sumpsf=sum(sum(sum(beadpsf,1),2),4);
Imidpsf=sumpsf(1,1,floor(length(zpos)/2)+1);
AA=2000; bg=40;
beadpsf=AA.*beadpsf./Imidpsf+bg;
beadpsfdx=AA.*beadpsfdx./Imidpsf+bg;
beadpsfdy=AA.*beadpsfdy./Imidpsf+bg;
beadpsfdz=AA.*beadpsfdz./Imidpsf+bg;

vIsum=beadpsf(:,:,floor(length(zpos)/2)+1,:)-bg;
vIsuma=sum(vIsum(:));
disp(vIsuma);


%% calculating the crlb for pixeled PSFs
if PSFsmooth==0
    ddxy=ddxy*1e3;
    ddz=ddz*1e3;
    du(:,:,:,:,1)=(beadpsfdx-beadpsf)./ddxy;
    du(:,:,:,:,2)=(beadpsfdy-beadpsf)./ddxy;
    du(:,:,:,:,3)=(beadpsfdz-beadpsf)./ddz;
    for kk=1:length(zpos)

        for ii=1:3
            for jj=1:3
                %tempM=(du(:,:,kk,ii).*du(:,:,kk,jj))./beadpsf(:,:,kk);
                tempM=zeros(size(du,1),size(du,2));
                for mm=1:2
                    tempM=tempM+1.*(du(:,:,kk,mm,ii).*du(:,:,kk,mm,jj))./(beadpsf(:,:,kk,mm)+1e-8);
                end
                tempM1=sum(sum(tempM,1),2);
                M(ii,jj)=tempM1(:);
            end
        end

        Minv=inv(M);
        CRLB(1:3,kk)=diag(Minv);

    end

    CRLB=sqrt(CRLB);
    figure(3)
    plot(CRLB(1,:),'r');
    hold on
    plot(CRLB(2,:),'r');
    title('lateral precisions by viPSF')
    figure
    plot(CRLB(3,:),'g');
    title('axial precison by viPSF')
end

