function imout = czt_getcztfft(fkxky,R,L,NN,MM)

%input
%fkxky: the pupil function
%R: half size of the pupil
%L: half size of the image
%NN: the number of pixels in pupil 1d
%MM: the number of pixels in image 1d

%output
% fourier transformed function
 
 %Bm=exp(-j*R*L*(1/NN-1)*(2/MM).*[0:MM-1]); Bm1 for x and Bm2 for y
 Bm=exp(-(j*2*pi*(1/NN-1)*2*L*(R/2/pi)/MM).*[0:MM-1]);
 Bm1=Bm.'*ones(1,NN);
 Bm2=Bm.'*ones(1,MM);  
 
 aczt=exp(j*R*L*(2/NN)*(1/MM-1));
 wczt=exp(-j*R*L*(4/(NN*MM)));
 %aczt=exp(j*2*pi*(1/MM-1)*2*L*(R/2/pi)/NN);
 %wczt=exp(-j*2*pi*4*L*(R/2/pi)/MM/NN);

 
 Etemp=(Bm1.*czt(fkxky,MM,wczt,aczt)).';
 imout=(Bm2.*czt(Etemp,MM,wczt,aczt)).';
 



end

