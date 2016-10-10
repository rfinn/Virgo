function iras60, lambda, flux

lambda0 = 60d0
readcol, 'z_response/iras60.txt', wave, resp, format='D,D',/silent

flux_wave = interpol(flux, lambda, wave)

fct_numerateur = flux_wave*resp/wave^2
int_numerateur = int_tabulated(wave, fct_numerateur)

fct_denumerateur = resp/(wave*lambda0)
int_denumerateur = int_tabulated(wave, fct_denumerateur)

result = int_numerateur/int_denumerateur

return, result
end
