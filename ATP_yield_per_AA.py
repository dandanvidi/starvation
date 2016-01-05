from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
import numpy as np
from scipy.integrate import quad


rho = 1100 # average cell density gr/liter
DW_fraction = 0.3 # fraction of DW of cells
Avogadro = 6.022 * 1e23 # Avogadro's number "exponent-less"
ATPM = 3.14 # mmol_gCDW_h
ATP_pool = 0.030 # mmol_gCDW from Bennett
glycogen_pool = (2.5/100)/0.18 # mmol_gCDW from Bennett
proteome_T = 10 # h proteome turnover

model = create_cobra_model_from_sbml_file('data/iJO1366.xml')
rxns = dict([(r.id, r) for r in model.reactions])
rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
rxns['ATPM'].lower_bound = 0
rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0
rxns['ATPS4rpp'].objective_coefficient = 1

AA = 'gly,ala_L,arg_L,asn_L,asp_L,cys_L,glu_L,gln_L,his_L,ile_L,leu_L,lys_L,met_L,\
phe_L,pro_L,ser_L,thr_L,trp_L,tyr_L,val_L'.split(',')

AA_to_symbol = {'gly':'G','ala_L':'A','arg_L':'R','asn_L':'N','asp_L':'D','cys_L':'C',
                'glu_L':'E','gln_L':'Q','his_L':'H','ile_L':'I','leu_L':'L','lys_L':'K',
                'met_L':'M','phe_L':'F','pro_L':'P','ser_L':'S','thr_L':'T','trp_L':'W',
                'tyr_L':'Y','val_L':'V'}
AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")
mmolATP_mmolAA = pd.DataFrame(index=AA_to_symbol.keys(), columns=['mol_ATP/mol_AA'])

aa_e = ''
for aa in AA:
    aa_e = 'EX_'+aa+'_e'
    rxns[aa_e].lower_bound = -1
    model.optimize()
    mmolATP_mmolAA.loc[aa] = model.solution.f
    rxns[aa_e].lower_bound = 0

mmolATP_mmolAA[mmolATP_mmolAA<1e-10] = 0
AA_mmol_gCDW = pd.DataFrame.from_csv('data/AA_mmol_gCDW.csv')['mmol/gCDW']
proteome_ATP_potential = mmolATP_mmolAA.mul(AA_mmol_gCDW, axis=0).sum().values[0] # mmol_gCDW

import matplotlib.pyplot as plt
time = np.arange(0,11)
proteome_decay = lambda a, x: (1-np.e**((-np.log(2)/proteome_T)*x))*a
missing = lambda x: x*ATPM - proteome_decay(proteome_ATP_potential, x)
energy_pool = glycogen_pool*25 + ATP_pool
t = 0
for i in np.arange(0,10,0.01):
    if quad(missing, 0, i)[0] > energy_pool:
        t = i
        break

limit = np.arange(0, t, 0.01)
plt.fill_between(limit, limit*ATPM, proteome_decay(proteome_ATP_potential, limit), color='0.7', 
                 edgecolor='', alpha=0.5, label='glycogen degradation')
plt.axvline(t, ls='-', color='k', lw=0.5, zorder=0)
plt.plot(time, proteome_decay(proteome_ATP_potential, time), label='proteome degradation')
plt.plot(time, time*ATPM, 'r', label='maintenance demand')
plt.xlabel('starvation duration [h]')
plt.ylabel('ATP [mmol/gCDW]')
plt.xticks(time)
#plt.yticks(np.arange(0,19,3))
#plt.ylim(0,30)
#plt.legend(loc=2, fontsize=10)
plt.tight_layout()
plt.savefig('ATP_from_proteome_degradation.png')