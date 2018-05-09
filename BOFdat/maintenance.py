"""
Maintenance
===========

Growth associated maintenance (GAM) is defined as the energy cost of growth, namely polymerization of macromolecules.
Non-Growth associated maintenance (NGAM) represent all the extra costs that the cell must overcome to operate.
This package offers two options for the user:
1- Calculate GAM and NGAM from experimental data
2- Estimate GAM and NGAM from theorical values

"""
def _import_model(path_to_model):
    import cobra
    extension = path_to_model.split('.')[-1]
    if extension == 'json':
        return cobra.io.load_json_model(path_to_model)
    elif extension == 'xml':
        return cobra.io.read_sbml_model(path_to_model)
    else:
        raise Exception('Model format not compatible, provide xml or json')

def _import_data(path_to_data):
    import pandas as pd
    import warnings

    data = pd.read_csv(path_to_data)
    '''
    #1- This file should have a header
    for i in data.columns:
        try:
            float(i)
            raise ValueError('Provide file header')
        except:
            pass
    '''
    return data

def experimental_maintenance(path_to_data, path_to_model,show_GAM=False):
    """

    Growth-associated maintenance (GAM) is the ATP cost of assembling macromolecules in the organism.
    This function calculates GAM from provided path to experimental data. This data includes growth rates on
    different carbon sources, the associated uptake rate for each carbon source and the secretion rates of metabolic
    wastes. More information on the format in which to provide the experimental data is available on GitHub.

    :param path_to_data: The data file is the outcome of the HPLC growth, uptake and secretion rate experiment.

    :param path_to_model: The path to the model, json or sbml formats supported

    :param show_GAM: bool, will associate colors with carbon sources for easier display later

    :return: a dictionary {GAM:value, NGAM:value}
    """

    #From experimental data growth rate on a given carbon source
    #Obtain best fit from model to experimental data
    def get_carbon_sources(data):
        return [c for c in data.Source]

    def attribute_colors(data,carbon_sources):
        # Attribute colors to carbon sources for ploting
        # Set a color palette
        import seaborn as sb
        color_palette = sb.color_palette('deep',len(carbon_sources))
        data['color'] = ''
        for i in len(carbon_sources):
            data.loc[data.Source == carbon_sources[i], 'color'] = color_palette[i]

        return data

    def get_biomass_objective_function(model):
        from cobra.util.solver import linear_reaction_coefficients
        return linear_reaction_coefficients(model).keys()[0]

    def atp_cost(model):
        solution = model.optimize().f
        if solution == None:
            solution = 0
        return solution

    def calculate_gam(model,data,show_GAM):
        import pandas as pd
        #Build the output matrix
        #Contains the ATP costs for each carbon source
        raw_GAM = pd.DataFrame(index=[data.loc[i,'Source'] for i in data.index],
                               columns=['Growth_rate', 'ATP', 'ATP_min', 'ATP_max', 'ATP_err']
                               )
        raw_GAM['Growth_rate'] = [data.loc[i,'GR'] for i in data.index]

        carbon_sources = get_carbon_sources(data)
        if show_GAM:
            data = attribute_colors(data,carbon_sources)
        else:
            data = data

        #Set parameters for all models
        #Set lower bound of BOF to 0
        #Get biomass objective function
        biomass = get_biomass_objective_function(model)
        biomass.lower_bound = 0.
        biomass.objective_coefficient = 0.

        # remove GAM from biomass function
        for key, value in biomass.metabolites.items():
            if abs(value) > 50:
                biomass.add_metabolites({key: -value})
        #Set carbon sources to 0
        '''
        Will need to prepare the model to run the simulations
        '''
        #model.reactions.EX_glc_LPAREN_e_RPAREN_.lower_bound = 0
        #Optimize for ATP maintenance
        model.reactions.ATPM.lower_bound = 0
        model.reactions.ATPM.objective_coefficient = 1.

        #model.reactions.EX_o2_LPAREN_e_RPAREN_.lower_bound = -1000
        atp, atp_min, atp_max, atp_err = [],[],[],[]
        for i in data.index:

            c_source = data.loc[i, 'Source']
            #Get growth rate
            growth_rate = data.loc[i, 'GR']
            growth_rate_err = data.loc[i, 'GR_std']
            #Get substrate uptake rate for the carbon source
            # gur = Glucose uptake rate
            substrate_uptake_rate = data.loc[i, 'SUR']
            # sur - substrate uptake rate
            substrate_uptake_rate_err = data.loc[i, 'SUR_std']

            # Mean
            #Excretion rates
            #Requires that the data file is generated appropriately
            for j in range(6,len(data.columns)):
                #Force flux through excretion
                if j % 2 == 0:
                    mean_excretion = data.loc[i,data.columns[j]]
                    model.reactions.get_by_id(data.columns[j]).lower_bound = mean_excretion

            '''
            Should return some errors if the file is not generated appropriately
            '''
            #Fix biomass at growth rate
            biomass.lower_bound = growth_rate

            # Iterate the Carbon source as %s
            model.reactions.get_by_id('%s' %(c_source,)).lower_bound = -1 * substrate_uptake_rate

            #Obtain solution and determine ATP cost
            atp.append(atp_cost(model))

            #Include experimental error in calculation
            # For both the excretion (max) and uptake (min) to get the max range
            #Minimum
            #Set biomass to growth rate
            biomass.lower_bound = growth_rate + growth_rate_err
            #Set uptake rate to the carbon source uptake rate
            model.reactions.get_by_id('%s' % c_source).lower_bound = -1 * (substrate_uptake_rate - substrate_uptake_rate_err)
            for j in range(6,len(data.columns)):
                #Force flux through excretion
                if j % 2 == 0:
                    mean_excretion = data.loc[i,data.columns[j]]
                    std_excretion = data.loc[i,data.columns[j+1]]
                    model.reactions.get_by_id(data.columns[j]).lower_bound = mean_excretion - std_excretion

            atp_min.append(atp_cost(model))

            #Maximum
            biomass.lower_bound = growth_rate - growth_rate_err
            for j in range(6,len(data.columns)):
                #Force flux through excretion
                if j % 2 == 0:
                    mean_excretion = data.loc[i,data.columns[j]]
                    std_excretion = data.loc[i,data.columns[j+1]]
                    model.reactions.get_by_id(data.columns[j]).lower_bound = mean_excretion + std_excretion

            atp_max.append(atp_cost(model))

            delta = atp_max[i] - atp_min[i]
            atp_err.append(delta / 2)

            #Reset uptake to 0
            model.reactions.get_by_id('%s' % c_source).lower_bound = 0.
            #Reset secretion to 0
            for j in range(6,len(data.columns)):
                #Force flux through excretion
                if j % 2 == 0:
                    mean_excretion = 0
                    model.reactions.get_by_id(data.columns[j]).lower_bound = mean_excretion

            '''
            Will have to define whether the aerobic/anaerobic option should be included
            '''
            # model.reactions.EX_o2_LPAREN_e_RPAREN_.lower_bound = -1000
        raw_GAM['ATP'] = atp
        raw_GAM['ATP_min'] = atp_min
        raw_GAM['ATP_max'] = atp_max
        raw_GAM['ATP_err'] = atp_err

        return raw_GAM

    def show_gam(raw_GAM):
        import seaborn as sns
        import numpy as np
        import matplotlib.pyplot as plt
        sns.set_style('whitegrid')
        x = raw_GAM['Growth_rate']
        y = raw_GAM['ATP']
        #Fit with np.polyfit
        m, b = np.polyfit(x, y, 1)

        print('m', m, 'b', b)

        #Get correlation
        x = raw_GAM['Growth_rate']
        y = raw_GAM['ATP']
        correlation = np.corrcoef(x, y)[0, 1]
        print('R2=', correlation ** 2)
        plt.scatter(raw_GAM['Growth_rate'], raw_GAM['ATP'])
        # plt.scatter(filtered_data['GR'],filtered_data['ATP'], color=filtered_data['color'], marker=filtered_data['marker'].tolist())
        plt.ylabel('ATP')
        plt.xlabel('Growth rate')
        plt.xlim([0, 1.1])
        plt.ylim([0, 110])
        plt.show()
        #plt.savefig('all_data.png')
        #plt.savefig('all_data.svg')
        plt.close()
        '''
        plt.errorbar(filtered_data['GR'], filtered_data['ATP'], xerr=filtered_data['GR_std'],
                     yerr=filtered_data['ATP_err'], fmt='.', ecolor='black')
        # plt.scatter(filtered_data['GR'],filtered_data['ATP'], color=filtered_data['color'])
        plt.plot(x, m * x + b, '-')

        plt.xlim([0, 1.1])
        plt.ylim([0, 120])
        plt.savefig('all_data_w_errors.png')
        plt.savefig('all_data_w_errors.svg')
        # embed()

        plt.show()
        plt.close()

        # plt.plot(x, y, '.')
        plt.scatter(filtered_data['GR'], filtered_data['ATP'], color=filtered_data['color'])
        plt.plot(x, m * x + b, '-')
        plt.xlim([0, 1.1])
        plt.ylim([0, 90])
        # plt.savefig('GAM_est.svg')
        plt.savefig('GAM_all_data_fit.png')
        plt.savefig('GAM_all_data_fit.svg')
        plt.show()
        plt.close()

        # for CI calc: combine dataset with std erors
        data2 = pd.DataFrame()
        data2 = data2.append(pd.DataFrame({'GR': filtered_data['GR'].tolist(), 'ATP': filtered_data['ATP'].tolist()}))

        # positive error
        test1 = filtered_data['GR'] + filtered_data['GR_std']
        test2 = filtered_data['ATP'] + filtered_data['ATP_err']
        data2 = data2.append(pd.DataFrame({'GR': test1, 'ATP': test2}))

        # negative error
        test1 = filtered_data['GR'] - filtered_data['GR_std']
        test2 = filtered_data['ATP'] - filtered_data['ATP_err']
        data2 = data2.append(pd.DataFrame({'GR': test1, 'ATP': test2}))

        # positive error
        # test1=filtered_data['GR']+filtered_data['GR_std']
        # test2=filtered_data['ATP']-filtered_data['ATP_err']
        # data2=data2.append(pd.DataFrame({'GR':test1,'ATP':test2}))

        # negative error
        # test1=filtered_data['GR']-filtered_data['GR_std']
        # test2=filtered_data['ATP']+filtered_data['ATP_err']
        # data2=data2.append(pd.DataFrame({'GR':test1,'ATP':test2}))



        sns.lmplot(x='GR', y='ATP', data=data2)
        plt.savefig('GAM_all_data_fit_CI.png')
        plt.savefig('GAM_all_data_fit_CI.svg')

        x = data2['GR']
        y = data2['ATP']

        # fit with np.polyfit
        m, b = np.polyfit(x, y, 1)
        print m, b

        plt.show()

        # embed()
        '''

        return {'GAM': m, 'NGAM': b}

    #1- Import model
    model = _import_model(path_to_model)
    #2- Import experimental data
    data = _import_data(path_to_data)
    #3- Calculate GAM
    raw_GAM = calculate_gam(model, data,show_GAM)
    #4-
    gams=show_gam(raw_GAM)
    #Grs = growth rates
    #new_gams_calc_grs(gams)

    return gams

def update_maintenance_costs(gams,model,RNA_atp):
    """
    This function updates the cost of maintenance in the BOF and the ATPM function.
    As of version 0.1.1 BOFdat assumes that the ATP maintenance reaction in the model is named ATPM and the
    biomass is set as the first objective of the model.

    :param gams: the dictionary generated by the experimental_maintenance function {GAM:x,NGAM:y}

    :param model: the model to update

    :param RNA_atp: the atp coefficient from the RNA

    """

    #Add the atp consumption to the biomass
    from BOFdat import update
    BOFdat.update.update_maintenance(gams,model,RNA_atp)

'''
Functions yet to be implemented
def theorical_GAM(protein_percent, rna_percent, dna_percent,CP = 4.324, CD = 1.365, CR = 0.406):
    #From experimental data: % of dry weight for each category
'''



