function phase_zern = zern_combine(zern_coef,zern_p)
    phase_zern=zeros(size(zern_p,1),size(zern_p,2));
    for ii=1:length(zern_coef)
        phase_zern=phase_zern+zern_coef(1,ii).*zern_p(:,:,ii);
    end
    phase_zern=exp(j.*phase_zern);
end