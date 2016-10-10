;-----------------------------------------------------------------------
function EEinv_dlumh, redshift
;-----------------------------------------------------------------------
common constantes_cosmo, hubble, Omega_matter, Omega_lambda, Omega_radiation

return, (Omega_matter*(1.+redshift)^3+Omega_lambda+Omega_radiation*(1.+redshift)^2)^(-0.5)

end

;-----------------------------------------------------------------------
function dlum_all, redshift,H0=H0, o_m=o_m, o_l=o_l
common constantes_cosmo, hubble, Omega_matter, Omega_lambda, Omega_radiation
;-----------------------------------------------------------------------

if not keyword_set(H0) then H0 = 70.
if not keyword_set(o_m) then o_m = 0.3
if not keyword_set(o_l) then o_l = 0.7

hubble=H0
c_light = 299792.48 ;km s^-1
DH= c_light/hubble ; Mpc
Omega_matter = o_m
Omega_radiation = 0.
Omega_lambda = o_l

result = (1.d0+redshift)*DH*qromb('EEinv_dlumh',0.,redshift)
return,result

end

