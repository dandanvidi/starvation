from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd

AA_data = pd.DataFrame.from_csv('data/AA_mmol_gCDW.csv')
prec = AA_data.precursor.str.split(' & ')
tmp = {aa:prec[aa] for aa in prec.index}
prec_sparse = {}
for k,v in tmp.iteritems():
    prec_sparse[k] = {}
    for m in v:
        m = m.split(' ')
        prec_sparse[k][m[1]] = int(m[0])
pd.DataFrame.from_dict(prec_sparse)
aa_to_prex = pd.DataFrame.from_dict(prec_sparse)
energy = ['ATP','NADH','NADPH']
prec_cost = pd.DataFrame.from_csv('data/precursors.csv')[energy]
aa_cost = AA_data[energy]
aa_cost_from_glucose = pd.DataFrame(index=AA_data.index, columns=energy)
for aa in AA_data.index:
    aa_cost_from_glucose.loc[aa] = prec_cost.mul(aa_to_prex[aa],axis=0).sum() + aa_cost.loc[aa]
aa_cost_from_glucose.to_csv('data/total_aa_cost_from_glucose.csv')



#model = create_cobra_model_from_sbml_file('data/iJO1366.xml')
#
#AA = 'gly,ala_L,arg_L,asn_L,asp_L,cys_L,glu_L,gln_L,his_L,ile_L,leu_L,lys_L,met_L,\
#phe_L,pro_L,ser_L,thr_L,trp_L,tyr_L,val_L'.split(',')
#
#mmolATP_mmolAA = pd.DataFrame(index=AA, columns=['cost'])
#rxns = {r.id:r for r in model.reactions}
#rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
#rxns['ATPM'].lower_bound = 0
#rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0
#rxns['ATPS4rpp'].objective_coefficient = 1
#aa_e = ''
#for aa in AA:
#    aa_e = 'EX_'+aa+'_e'
#    rxns[aa_e].lower_bound = -1
#    model.optimize()
#    mmolATP_mmolAA.loc[aa] = model.solution.f
#    rxns[aa_e].lower_bound = 0
#
#mmolATP_mmolAA[mmolATP_mmolAA<1e-10] = 0
#
#for aa in AA: