function [kx,xx] = czt_getcztcoord(num_x,NN1)


%input 
% num_x, the number of pixels in the image plane
%NN1,minimum pixels in the kx
%output
%kx xx, between [-1, 1] or [-1, 1]

if nargin<1
NN1=500;
end

MM=num_x;

% sample for dk|z|<pi 【-N/2，N/2】代表k取值不为0的区域
%NN2=6*(num_z*dz/2)*2*na^2/(sqrt(n0^2-na^2));
NN=max([NN1,num_x]);
NN=2^(floor(log2(NN+MM))+1)-MM+1;
%dk=2*kmax/NN;     % 
 

% R=dk*NN/2 ;   u0=dx*MM/2;
 kx=linspace(-1+1/NN,1-1/NN,NN);   
 xx=linspace(-1+1/MM,1-1/MM,MM);
 

end

