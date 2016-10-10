function chary_elbaz_24um,redshift,SmicroJy, save=save, Ho=Ho, o_m=o_m, o_l=o_l,verbose=verbose

if not keyword_set(Ho) then Ho = 70.
if not keyword_set(o_m) then o_m = 0.3
if not keyword_set(lambda) then lambda = 0.7

ngals=n_elements(redshift) 
if ngals ne n_elements(SmicroJy) then begin
    print,'ERROR in dimensions'
    return,0
endif

distance=dlum_all(redshift, H0=H0, o_m=o_m, o_l=o_l)
;
lambda_0 = 23.675               ; microns
c_light_um = 2.99793d14         ; microns/s
;
restore,'chary_elbaz.save'
;
Slambda_seds_obs=dblarr(105,1366)
Snu_seds_obs=dblarr(105,1366)
;
readcol, 'z_response/mips24lg.dat', l_r24, r_r24, format='D,D', /silent
NN_r24 = n_elements(l_r24)
Slambda_seds_24um = dblarr(105, NN_r24)

T_corps_noir = 10000            ; K
corps_noir = planck(l_r24*1d4, T_corps_noir)
num_seds = dblarr(105)
denom_seds = int_tabulated(l_r24, r_r24*corps_noir)
Slambda_seds_mips24 = dblarr(105)


Snu = dblarr(ngals,1366)
luminosity =  dblarr(ngals,1366)
LIR_Sanders=dblarr(ngals)
SFR_Sanders=dblarr(ngals)
S24um=dblarr(ngals)

for igal=0,ngals-1 do begin

    if (keyword_set(verbose) and igal mod 10 eq 0) then print,igal,'/',ngals

    ;; SmicroJy in 10^-29 erg/s/cm^2/Hz
    Slambda = 1d-29*1d-4*(c_light_um/lambda_0^2)*SmicroJy(igal) ; in erg/s/cm^2/Angstrom

    lambda_obs   = lambda*(1d0+redshift(igal))

    for i=0,104 do Snu_seds_obs(i,*) = (1d0+redshift(igal))*nuLnuinLsun(i,*)/(1d-32*4*!dpi*(distance(igal)*3.0856d22)^2*(3d14/lambda)/3.826d26) ; in microJy

    for i=0,104 do Slambda_seds_obs(i, *)= 3.826d33*nuLnuinLsun(i,*)/(4*!dpi*(distance(igal)*3.0856d24)^2*(1d4*lambda))/(1d0+redshift(igal)) ; erg/s/cm^2/Angstrom
;
    ;; spectrum of  each template interpolated at the 24um response curve wavlelengths
    for i=0,104 do Slambda_seds_24um(i, *) = interpol(Slambda_seds_obs(i, *), lambda_obs, l_r24) 

    ;; integration over the 24 response curve
    for i=0,104 do num_seds(i) = int_tabulated(l_r24, r_r24*Slambda_seds_24um(i, *))

    ;; to scale to a black body
    norm_seds = num_seds/denom_seds 

    ;; model flux density at 24um (over the filter)
    for i=0,104 do Slambda_seds_mips24(i) = norm_seds(i)*interpol(corps_noir, l_r24*1d4, lambda_0*1d4) ; in erg/s/cm^2/Angstrom

    ind_sed=interpol(findgen(105),Slambda_seds_mips24,Slambda) ;; get index for this observed flux
    if (ind_sed lt 0 or ind_sed gt 104) then begin
        print,'Warning: attempt to extrapolate the 24um luminosity outside the model limits !!!'
    endif

    Snu(igal,*)=reform(interpolate(Snu_seds_obs,ind_sed,findgen(1366),/grid))

    f12_rf = iras12(lambda, Snu(igal,*)/(1d0+redshift(igal)))
    f25_rf = iras25(lambda, Snu(igal,*)/(1d0+redshift(igal)))
    f60_rf = iras60(lambda, Snu(igal,*)/(1d0+redshift(igal)))
    f100_rf = iras100(lambda, Snu(igal,*)/(1d0+redshift(igal)))


    Lir_sanders(igal) = 1.8d-14*1d-6*(13.48*f12_rf + 5.16*f25_rf + 2.58*f60_rf + f100_rf)*4*!dpi*(distance(igal)*3.0856d22)^2/3.826d26
    SFR_sanders(igal) = 1.7217d-10*Lir_sanders(igal)

    S24um(igal) = SmicroJy(igal)
endfor

if keyword_set(save) then save, filename='chary_elbaz_24um.save', Ho, redshift, S24um, Lir_sanders, SFR_sanders,  lambda, luminosity, Snu

return, Lir_sanders

end
