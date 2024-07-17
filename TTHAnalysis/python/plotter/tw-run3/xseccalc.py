# -*- coding: utf-8 -*-

# ========================================================
# Datos de entrada
mu     = 0.936
muup   = 0.121
mudn   = 0.118

mustat   = 0.024
musystup = 0.113
musystdn = 0.110
mulumi   = 0.038

sigmath = 87.9

ndectot = 1 #### AJUSTAR ESTO EN FUNCIÓN A LAS CIFRAS SIGNIFICATIVAS Y A LAS GUIDELINES DE CMS DE PUBLICACIÓN
ndecsep = 1 #### AJUSTAR ESTO EN FUNCIÓN A LAS CIFRAS SIGNIFICATIVAS Y A LAS GUIDELINES DE CMS DE PUBLICACIÓN

# ========================================================
# Calculos
totup = (mustat**2 + musystup**2 + mulumi**2)**(1./2)
totdn = (mustat**2 + musystdn**2 + mulumi**2)**(1./2)

muup_rel = muup/mu
mudn_rel = mudn/mu

mustat_rel   = mustat/mu
musystup_rel = musystup/mu
musystdn_rel = musystdn/mu
mulumi_rel   = mulumi/mu

sigma   = sigmath * mu
sigmaup = muup_rel*sigma
sigmadn = mudn_rel*sigma
sigmastat   = mustat_rel*sigma
sigmasystup = musystup_rel*sigma
sigmasystdn = musystdn_rel*sigma
sigmalumi   = mulumi_rel*sigma

sigma_tot   = round(sigma,   ndectot)
sigmaup_tot = round(sigmaup, ndectot)
sigmadn_tot = round(sigmadn, ndectot)

sigma_sep       = round(sigma,       ndecsep)
sigmastat_sep   = round(sigmastat,   ndecsep)
sigmasystup_sep = round(sigmasystup, ndecsep)
sigmasystdn_sep = round(sigmasystdn, ndecsep)
sigmalumi_sep   = round(sigmalumi,   ndecsep)
sigmaup_sep     = (sigmastat_sep**2 + sigmasystup_sep**2 + sigmalumi_sep**2)**(1./2)
sigmadn_sep     = (sigmastat_sep**2 + sigmasystdn_sep**2 + sigmalumi_sep**2)**(1./2)


# ========================================================
# Comprobacion incs separadas y juntas
print("\n COMPROBACIÓN: total vs. suma cuadrática de stat, lumi y syst (en absoluto)")
print("Total arriba: \t" + str(muup) + ", total sumado arriba: \t" + str(totup))
print("Total abajo: \t"  + str(mudn) + ", total sumado abajo: \t"  + str(totdn))


print("\n Fuerza de señal - Resultados absolutos")
print("mu = " + str(mu) + " +\t" + str(muup) + "\t- " + str(mudn))
print("mu = " + str(mu) + " +-\t" + str(mustat) + " (stat.)\t+ " + str(musystup) + " - " + str(musystdn) + " (syst.) +- " + str(mulumi) + " (lumi.)")

print("\n Fuerza de señal - Resultados relativos")
print("mu = " + str(mu) + " +\t" + str(round(muup_rel*100, 4)) + "% \t- " + str(round(mudn_rel*100, 4)) + " %")
print("mu = " + str(mu) + " +-\t" + str(round(mustat_rel*100, 4)) + "% (stat.)\t+ " + str(round(musystup_rel*100, 4)) + " - " + str(round(musystdn_rel*100, 4)) + "% (syst.) +- " + str(round(mulumi_rel*100, 4)) + "% (lumi.)")

print("\n Sección eficaz - Resultados absolutos")
print("sigma = " + str(round(sigma, 4)) + " +  "  + str(round(sigmaup, 4)) + " - " + str(round(sigmadn, 4)) + " pb")
print("sigma = " + str(round(sigma, 4)) + " +- " + str(round(sigmastat, 4)) + " (stat.)\t+ " + str(round(sigmasystup, 4)) + " - " + str(round(sigmasystdn, 4)) + " (syst.) +- " + str(round(sigmalumi, 4)) + " (lumi.) pb")

print("\n COMPROBACIÓN: resultados finales redondeando a cifras significativas siguiendo la guía de estilo de CMS, suma cuadrática, y su relativo")
print("sigma = " + str(sigma_tot) + " +  "  + str(sigmaup_tot) + " - " + str(sigmadn_tot) + " pb; + " + str(sigmaup_tot/sigma_tot*100) + " - " + str(sigmadn_tot/sigma_tot*100) + " %") 
#print "sigma = " + str(sigma_sep) + " +- " + str(sigmastat_sep) + " (stat.)\t+ " + str(sigmasystup_sep) + " - " + str(sigmasystdn_sep) + " (syst.) +- " + str(sigmalumi_sep) + " (lumi.) pb; + " + str(sigmaup_sep) + " - " + str(sigmadn_sep) + " pb; + " + str(sigmaup_sep/sigma_sep*100) + " - " + str(sigmadn_sep/sigma_sep*100) + " %"
print("sigma = " + str(sigma_sep) + " +- " + str(sigmastat_sep) + " (stat.)\t+ " + str(sigmasystup_sep) + " - " + str(sigmasystdn_sep) + " (syst.) +- " + str(sigmalumi_sep) + " (lumi.) pb; + " + str(sigmaup_sep/sigma_sep*100) + " - " + str(sigmadn_sep/sigma_sep*100) + " %")

