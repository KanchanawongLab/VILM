function crlbout= obj_ipalmpsfoptim(fit_par_short)


global showPSFs  % set to 0, PSFs and CRLB will be displayed.
global addvortex;
global kxx1 kyy1 kz_c1 kz_c2
global CTF1
global R L NN MM 
global zpos
global zern_p

global k0;
global AA bg

sz=length(zpos);
mid_sz=floor((sz+1)/2);
zern_coef=fit_par_short(1:end);  
phase_zern=zern_combine(zern_coef,zern_p);  % generate zernike phase 

% initial parameter for calculating psfs
x0=0;y0=0;z=0;
%focus shift when considering the top objective
z1=0;z2=0;

%parameters for calculating crlb
ddxy=0.01; % 10nm
ddz=1e-3;  %1nm

%whether to add vortex phase mask
fai=atan2(kxx1,kyy1);
if addvortex==1  
    fmask2=exp(j.*fai);
    fmask3=exp(-j.*fai);
else
    fmask2=0.*phase_zern+1;
    fmask3=fmask2;
end

%% CALCULATE THE PSFs
for jj=1:sz  
    focpos=zpos(jj);
   
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z,focpos);
    fmask1=phase_zern;
    fkxky=fpupil.*fmask1.*fmask2;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z1,z2-focpos);
    fmask1=exp(-j.*angle(phase_zern));
    fkxky=fpupil.*fmask1.*fmask3;
    im1=czt_getcztfft(fkxky,R,L,NN,MM); 
    beadpsf(:,:,jj,1)=abs(im+j.*im1).^2;
    beadpsf(:,:,jj,2)=abs(im-j.*im1).^2;
    
    % figure(12)
    % imagesc(imag(phase_zern));
    %
    if showPSFs==1
        h1=figure(1);
        set(h1,'Position',[100 100 600 300])
        subplot(121)
        imagesc(abs(im+j.*im1).^2);
        title('PSF in channel one')
        colormap('copper')
        subplot(122)
        imagesc(abs(im-j.*im1).^2);
        title('PSF in channel two')
        colormap('copper')
    end
    
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+ddxy,y0,z,focpos);
    fmask1=phase_zern;
    fkxky=fpupil.*fmask1.*fmask2;
    im=czt_getcztfft(fkxky,R,L,NN,MM);  
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0+ddxy,y0,z,-focpos);
    fmask1=exp(-j.*angle(phase_zern));
    fkxky=fpupil.*fmask1.*fmask3;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);   
    beadpsfdx(:,:,jj,1)=abs(im+j.*im1).^2;
    beadpsfdx(:,:,jj,2)=abs(im-j.*im1).^2;
    
    
    
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0+ddxy,z,focpos);
    fmask1=phase_zern;
    fkxky=fpupil.*fmask1.*fmask2;
    im=czt_getcztfft(fkxky,R,L,NN,MM); 
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0+ddxy,z,-focpos);
    fmask1=exp(-j.*angle(phase_zern));
    fkxky=fpupil.*fmask1.*fmask3;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsfdy(:,:,jj,1)=abs(im+j.*im1).^2;
    beadpsfdy(:,:,jj,2)=abs(im-j.*im1).^2;
    
    
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z,focpos+ddz);
    fmask1=phase_zern;
    fkxky=fpupil.*fmask1.*fmask2;
    im=czt_getcztfft(fkxky,R,L,NN,MM);
    fpupil=CTF1.*czt_getphasexyz(kxx1,kyy1,kz_c1,kz_c2,x0,y0,z,-focpos-ddz);
    fmask1=exp(-j.*angle(phase_zern));
    fkxky=fpupil.*fmask1.*fmask3;
    im1=czt_getcztfft(fkxky,R,L,NN,MM);
    beadpsfdz(:,:,jj,1)=abs(im+j.*im1).^2;
    beadpsfdz(:,:,jj,2)=abs(im-j.*im1).^2;
    
    
end


%% normalize the psf, multiply with photon number, and add background noise
sumpsf=sum(sum(sum(beadpsf,1),2),4);
Imidpsf=sumpsf(1,1,mid_sz);
AA=2000; bg=40;
beadpsf=AA.*beadpsf./Imidpsf+bg;
beadpsfdx=AA.*beadpsfdx./Imidpsf+bg;
beadpsfdy=AA.*beadpsfdy./Imidpsf+bg;
beadpsfdz=AA.*beadpsfdz./Imidpsf+bg;


%% calculating crlb
crlb_out=0;
ddxy=1e3*ddxy; ddz=1e3*ddz;
du(:,:,:,:,1)=(beadpsfdx-beadpsf)./ddxy;
du(:,:,:,:,2)=(beadpsfdy-beadpsf)./ddxy;
du(:,:,:,:,3)=(beadpsfdz-beadpsf)./ddz;
for kk=1:sz 
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

%% calculate the weighted CRLB
crlbout=sum(CRLB(:));
crlbout=crlbout/120+max(CRLB(3,:));

%% monitor whether it decreased during the optimization
disp(crlbout)

if showPSFs==1
    figure(3)
    plot(CRLB(1,:));
    hold on
    plot(CRLB(2,:))
    plot(CRLB(3,:))
    hold off
end

end
