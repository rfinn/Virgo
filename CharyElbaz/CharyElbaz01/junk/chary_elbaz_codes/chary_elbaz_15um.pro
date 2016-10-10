function chary_elbaz_15um,redshift,SmicroJy, save=save, Ho=Ho, o_m=o_m, o_l=o_l

if not keyword_set(Ho) then Ho = 70.
if not keyword_set(o_m) then o_m = 0.3
if not keyword_set(o_l) then o_l = 0.7

distance=dlum_all(redshift, H0=H0, o_m=o_m, o_l=o_l)
;
readcol,'z_response/lw3.response',rl,rr,format='D,D',/silent
norm=int_tabulated(rl,rr)
;
restore,'chary_elbaz.save'
;
Snu_seds_obs=dblarr(105,1366)
;
lambda_obs   = lambda*(1d0+redshift)
for i=0,104 do Snu_seds_obs(i,*) = (1d0+redshift)*nuLnuinLsun(i,*)/(1d-32*4*!dpi*(distance*3.0856d22)^2*(3d14/lambda)/3.826d26) ; microJy
;
Snu_obs = dblarr(105,n_elements(rl))
for i=0,104 do Snu_obs(i,*) = interpol(Snu_seds_obs(i,*),lambda_obs,rl)
;
Snu_seds_lw3 = dblarr(105)
for i=0,104 do Snu_seds_lw3(i) = int_tabulated(rl,Snu_obs(i,*)*14.3*rr/rl)/norm
;
; 14.3 microns <-> lambda(262)
;
choice_sed = min(abs(SmicroJy - Snu_seds_lw3),ind_sed)
corr = SmicroJy/Snu_seds_lw3(ind_sed)

luminosity_interm = corr*nuLnuinLsun(ind_sed,*)
luminosity =  dblarr(1366)
luminosity(*) =  luminosity_interm(0, *)
NN = 100
lambda_ir = 8d0*(1d3/8d0)^(indgen(NN)/(NN-1d0))
nuLnu_ir  = interpol(luminosity,lambda,lambda_ir)
Lir = int_tabulated(lambda_ir,nuLnu_ir/lambda_ir)
SFR = 1.7217d-10*Lir

Snu = dblarr(1366)
Snu(*) =  Snu_seds_obs(ind_sed, *)*corr

f12 = iras12(lambda, Snu/(1d0+redshift))
f25 = iras25(lambda, Snu/(1d0+redshift))
f60 = iras60(lambda, Snu/(1d0+redshift))
f100 = iras100(lambda, Snu/(1d0+redshift))
Lir_sanders = 1.8d-14*1d-6*(13.48*f12 + 5.16*f25 + 2.58*f60 + f100)*4*!dpi*(distance*3.0856d22)^2/3.826d26
SFR_sanders = 1.7217d-10*Lir_sanders

return, Lir_sanders
end
