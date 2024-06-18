function phase_xyz=czt_getphasexyz(kxx,kyy,w1,w2,x0,y0,z,focpos)

%add the transverse shift term
phase_xy=exp(j.*x0.*kxx+j.*y0.*kyy);

% add the axial shift term
phase_xyz=phase_xy.*exp(j.*real(w2).*z-j.*real(w1).*focpos).*exp(-imag(w2)*z);




end

