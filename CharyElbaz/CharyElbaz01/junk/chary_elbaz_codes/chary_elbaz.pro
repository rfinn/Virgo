function chary_elbaz,redshift,SmicroJy, lam_filter, save=save, H0=H0, o_m=o_m, o_l=o_l, resol=resol

if not keyword_set(H0) then H0 = 70.
if not keyword_set(o_m) then o_m = 0.3
if not keyword_set(lambda) then lambda = 0.7
if not keyword_set(resol) then Resol=3.

distance=dlum_all(redshift, H0=H0, o_m=o_m, o_l=o_l)

restore,'chary_elbaz.save'
;
Snu_seds_obs=dblarr(105,1366)
;
lambda_obs   = lambda*(1d0+redshift)
for i=0,104 do Snu_seds_obs(i,*) = (1d0+redshift)*nuLnuinLsun(i,*)/(1d-32*4*!dpi*(distance*3.0856d22)^2*(3d14/lambda)/3.826d26) ; microJy
;
;
Snu_seds_instrument = dblarr(105)
;
lam_min=lam_filter-lam_filter/(2.*Resol)
lam_max=lam_filter+lam_filter/(2.*Resol)

for i=0,104 do Snu_seds_instrument(i) = mean(Snu_seds_obs(i,where(lambda_obs ge lam_min and lambda_obs le lam_max)))

choice_sed = min(abs(SmicroJy - Snu_seds_instrument),ind_sed)

Snu = dblarr(1366)
Snu(*) =  Snu_seds_obs(ind_sed, *)*(SmicroJy/Snu_seds_instrument(ind_sed))

f12 = iras12(lambda, Snu/(1d0+redshift))
f25 = iras25(lambda, Snu/(1d0+redshift))
f60 = iras60(lambda, Snu/(1d0+redshift))
f100 = iras100(lambda, Snu/(1d0+redshift))
Lir_sanders = 1.8d-14*1d-6*(13.48*f12 + 5.16*f25 + 2.58*f60 + f100)*4*!dpi*(distance*3.0856d22)^2/3.826d26

return, Lir_sanders

end
