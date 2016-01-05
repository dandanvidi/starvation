# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:54:11 2016

@author: dan
"""

from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
from cobra.core.Formula import Formula

model = create_cobra_model_from_sbml_file('data/iJO1366.xml')
rxns = dict([(r.id, r) for r in model.reactions])
rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized 

exchange_reactions = [r for r in model.reactions if 'EX_' in r.id]
carbon_sources = pd.DataFrame(index=model.reactions, 
                              columns=['metabolite id', 'metabolite name', 'MW [Da]', 
                                       'carbons', 'yield [gCDW/mmol carbon]',
                                       'ATP yield [mol ATP/mol carbon]', 'exchange'])
                                       
for i, r in enumerate(exchange_reactions):
    print i, 
    m = r.metabolites.keys()[0]
    mf = Formula(m.formula)
    if r.lower_bound == 0 and 'C' in mf.elements: # potential carbon source
        uptake = 1. / mf.elements['C']
        r.lower_bound = -uptake
        model.optimize()
        
        if model.solution.f: # E. coli can grow on this carbon source
            print m.name
            fluxes = model.solution.x_dict

            carbon_sources['metabolite id'][r] = m.id
            carbon_sources['metabolite name'][r] = m.name
            carbon_sources['MW [Da]'][r] = mf.weight            
            carbon_sources['carbons'][r] = mf.elements['C']            
            carbon_sources['yield [gCDW/mmol carbon]'][r] = model.solution.f
                                             
            ex = {}
            for r0 in exchange_reactions:
                m0 = r0.metabolites.keys()[0]
                m0f = Formula(m0.formula)
                if 'C' in m0f.elements > 0 and fluxes[r0.id] > 0:
                    ex[r0.id] = fluxes[r0.id]
            
            carbon_sources['exchange'][r] = str(ex)
            
            # calculate ATP yiled
            rxns['ATPM'].lower_bound = 0
            rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0
            rxns['ATPS4rpp'].objective_coefficient = 1
            model.optimize()

            rxns['ATPS4rpp'].objective_coefficient = 0
            rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 1

            carbon_sources['ATP yield [mol ATP/mol carbon]'][r] = model.solution.f
                  
        r.lower_bound = 0            
        
               
carbon_sources.dropna(how='all', inplace=True)
carbon_sources.to_csv('res/carbon_exchange_reactons.tsv', sep='\t')

#        carbon_sources[r.id] = r.metabolites.keys()[0]
#

#rxns['ATPM'].lower_bound = 0
#rxns['Ec_biomass_iJO1366_core_53p95M'].objective_coefficient = 0
#rxns['ATPS4rpp'].objective_coefficient = 1