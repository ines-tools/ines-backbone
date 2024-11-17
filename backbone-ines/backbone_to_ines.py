import os
import sys
import yaml
import numpy as np
from datetime import datetime
from pathlib import Path
import math
#from ines_tools import ines_transform
import ines_transform
import spinedb_api as api
from sqlalchemy.exc import DBAPIError
from spinedb_api.exception import NothingToCommit
#from spinedb_api import purge

def main():
    # transform spine db with backbone data (source db) into a spine db that already has the ines structure (target_db)
    with api.DatabaseMapping(url_db_in) as source_db:
        with api.DatabaseMapping(url_db_out) as target_db:
            # completely empty database
            #purge.purge(target_db, purge_settings=None)
            # add ines structure
            # empty database except for ines structure
            target_db.purge_items('parameter_value')
            target_db.purge_items('entity')
            target_db.purge_items('alternative')
            target_db.refresh_session()
            target_db.commit_session("Purged everything except for the existing ines structure")
            # copy alternatives and scenarios
            for alternative in source_db.get_alternative_items():
                target_db.add_alternative_item(name=alternative["name"])
            for scenario in source_db.get_scenario_items():
                target_db.add_scenario_item(name=scenario["name"])
            for scenario_alternative in source_db.get_scenario_alternative_items():
                target_db.add_scenario_alternative_item(
                    alternative_name=scenario_alternative["alternative_name"],
                    scenario_name=scenario_alternative["scenario_name"],
                    rank=scenario_alternative["rank"]
                )
            # commit changes
            target_db.refresh_session()
            target_db.commit_session("Added scenarios and alternatives")
            #Add entities that are always there
            target_db = add_base_entities(target_db)

            target_db, timestep_names = create_timeline(source_db, target_db)

            # copy entities from yaml files
            print("copy entities")
            target_db = ines_transform.copy_entities(source_db, target_db, entities_to_copy)
            #create in-and-out relationships
            print("add link capacities")
            target_db = process_links(source_db, target_db)
            print("create unit relationships")
            target_db = create_unit_relationship(source_db, target_db)
            #copy parameters to relationships
            print("create relationships from parameters")
            target_db = ines_transform.transform_parameters_to_relationship_entities(source_db, target_db, parameters_to_relationships)
            #copy parameters to entities, but the entity name is from a parameter
            print("add direct parameters to different name")
            target_db = ines_transform.transform_parameters_entity_from_parameter(source_db, target_db, parameters_to_parameters)   
            # copy numeric parameters
            print("add direct parameters")
            target_db = ines_transform.transform_parameters(source_db, target_db, parameter_transforms)
            # copy method parameters
            print("add process methods")
            target_db = ines_transform.process_methods(source_db, target_db, parameter_methods)
            print("add capacities")
            target_db = process_capacities(source_db, target_db)
            print("create timeseries from price change")
            target_db = create_price_change(source_db, target_db, timestep_names)
            print("create profiles")
            target_db = create_profiles(source_db, target_db)
            print("diffusion coefficient")
            target_db = diff_coeff(source_db, target_db)
            print("add node types")
            target_db = add_node_types(source_db, target_db)
            # copy entities to parameters
            #target_db = ines_transform.copy_entities_to_parameters(source_db, target_db, entities_to_parameters)

            # manual scripts
            # copy capacity specific parameters (manual scripting)
            #target_db = process_capacities(source_db, target_db)
            try:
                target_db.commit_session("Added entities")
            except NothingToCommit:
                pass
            except DBAPIError as e:
                print("failed to commit entities and entity_alternatives")
            return target_db

# only the part below is specific to a tool

def get_settings():
    convertpath = 'backbone_to_ines_settings.yaml'
    with open(convertpath,'r') as file:
        settings = yaml.safe_load(file)

    return settings

# quick conversions using dictionaries
# these definitions can be saved here or in a yaml configuration file
'''
    conversion_configuration

A function that saves/loads from yaml files (currently only supported file type). The data is also available within this function but is only loaded when requested.

If a filepath is given and it exists, it will be loaded. If it does not exist, data from within this function will be saved to the file (if available).

If a filename is given, the data from this function will be returned.

conversions : list of file paths or file names
overwrite : boolean that determines whether an existing file is overwritten with the data inside this function

return a list of conversion dictionaries
'''
def conversion_configuration(conversions = ['backbone_to_ines_entities', 'backbone_to_ines_parameters','backbone_to_ines_parameter_methods',
                                            'backbone_to_ines_parameters_to_relationships','backbone_to_ines_parameters_to_parameters'], overwrite=False):
    returnlist = []
    for conversion in conversions:
        # default is data from within this function
        convertname = conversion
        load = False
        save = False

        # check whether a file or name is passed and reconfigure this function accordingly
        convertpath = Path(conversion)
        if convertpath.suffix == '.yaml':
            convertname = convertpath.stem
            if convertpath.is_file() and not overwrite:
                load = True
            else:
                save = True

        if load:
            # load data from file
            with open(convertpath,'r') as file:
                returnlist.append(yaml.safe_load(file))
        else:
            # get data from within this function
            convertdict = None
            if convertname == 'backbone_to_ines_entities':
                convertdict = {
                    #boundary = ?
                    'constraint': ['constraint'],
                    #grid: group?
                    'group': ['set'],
                    'restype': ['reserve'],
                    'unit': ['unit'],
                    #up_down: add to ines?
                    'node': ['node'],
                }
            if convertname == 'backbone_to_ines_parameters':
                convertdict = {
                    'constraint': {
                        #'constraint': {
                        #},
                    },
                    'group':{
                        'set':{
                            #ROCOF:  #inertia stuff
                            'constrainedCapTotalMax': 'invest_max_total', #check units
                            #'constrainedOnlineTotalMax':   not in ines
                            #'defaultFrequency':  inertia stuff
                            #'dynamicInertia' # not in ines
                            'energyMax': 'flow_max_cumulative',
                            'energyMin':'flow_min_cumulative', 
                            #'energyShareMax'
                            #'energyShareMin'
                            #'staticInertia' part of the ROCOF stuff
                        }
                    },
                    'model':{ 
                        #'system':{
                            #'discountFactor': not 'discount_rate'! calculate discount rate? period wise discount rate??
                            #'t_invest': #assumed to be the period start in ines
                        #}
                    },
                    'node':{
                        'node':{
                            'price': 'commodity_price'
                        }
                    },
                    'restype':{
                        #'restype':{
                            #'restypeReleasedForRealization'
                            #'restype_inertia'
                        #}
                    },
                    'unit':{
                        'unit':{
                            'availability': 'availability',
                            #'becomeAvailable':  not ines or calculate to the availability timeseries
                            #'becomeUnavailable': not ines or calculate to the availability timeseries
                            #'boundSamples': not in ines
                            #'efficiency':  #eff curve calculate
                            #'efficiency_ts':  #eff curve calculate
                            #'initialOnlineStatus': #hot start not in ines
                            #'investMIP': #method for MIP invest
                            'maxUnitCount': 'units_max_cumulative',
                            'minOperationHours': ['min_uptime', 60],
                            'minShutDownHours': ['min_downtime', 60],
                            'minUnitCount': 'units_min_cumulative',
                            #'rampSpeedFromMinLoad': not in ines
                            #'rampSpeedToMinLoad' not in ines
                            #'startColdAfterXhours'
                            #'startWarmAfterXHours'
                            'unitCount': 'units_existing',
                            #'unit_fail': #reserve method calculation
                            #'useInitialOnlineStatus' #ines hot start flag
                        }
                    },
                    'grid__node':{
                        'node':{
                            #'capacityMargin' #not in ines
                            #'energyStoredPerUnitOfState': #coefficient to the state params
                            'influx': ['flow_profile',1.0,[[2]]]
                            #'selfDischargeLoss: not in ines
                            #'storageValue': 'storagePrice' not ines
                            #'storageValueUseTimeSeries' flag
                            #'usePrice' flag for priceChange
                        }
                    },
                    'group__emission':{
                        #'set': {
                            #'emissionCap': ['co2Total', 1.0, [1]] #might work
                            #'emissionPriceChange': #create timeseries?
                            # 'emissionTax': ['co2Price', 1.0, [1]] #might work
                        #}
                    },
                    'group__restype':{
                        #gate_closure
                        #reserve_activation_duration
                        #reserve_length
                        #reserve_reactivation_time
                        #update_frequency
                        #update_offset
                        #useTimeseries flag
                    },
                    'node__emission':{
                        #'node':{
                            #'emission_content': 'co2_content' # check emission, add SO2
                        #}
                    },
                    'unit__group':{
                        #'constrainedCapMultiplier': #for constrainedCapTotalMax not in ines 
                        #'constrainedOnlineMultiplier':  #for constrainedOnlineTotalMax not in ines 
                    },
                    'unit__node':{
                        #'fixedFuelFraction': #Fixed share of a fuel in the start-up fuel consumption mix
                    },
                    'grid__node__boundary':{
                        #constant
                        #multiplier
                        #slackCost
                        #timeseries
                        #useConstant
                        #useTimeseries
                    },
                    'grid__node__group':{
                    },
                    'grid__restype__up_down':{
                        #'LossOfTrans': n-1 flag
                        #'reserveDemand':
                    },
                    'grid__node__unit__boundary':{
                        #'rampCost'
                        #'rampLimit'
                    },
                    'grid__node__unit__emission':{
                        #'fomEmissions': fixed emissions  
                        #'invEmissionFactor':   investment emissions possibility to divide over years
                        #'invEmissions':    emissions from investments
                        #'vomEmissions':    variable emissions
                    },
                    'grid__node__unit__restype':
                        {
                            #down
                            #offlineReserveCapability
                            #portion_of_infeed_to_reserve
                            #reserveReliability
                            #reserve_increase_ratio
                            #up
                        }
                    ,'grid__node__unit__io':{
                        'unit':{
                            'shutdownCost': ['shutdown_cost',1.0,[[2]]]
                        }
                    }    
                }
            if convertname == 'backbone_to_ines_parameter_methods':
                convertdict = {
                    'grid__node': {
                        'node':{
                            'boundAll':{
                                1.0:{
                                    'storage_state_fix_method': ['fix_all_timesteps', [2]]
                                }
                            }, 
                            'boundEnd':{
                                1.0:{
                                    'storage_state_fix_method': ['fix_horizon_end', [2]]
                                }
                            },
                            'boundStart':{
                                1.0:{
                                    'storage_state_fix_method': ['fix_start_and_horizon_end', [2]]
                                }
                            },
                            'boundStartAndEnd':{
                                1.0:{
                                    'storage_state_fix_method': ['fix_start', [2]]
                                }
                            },
                            'boundSamples':{
                                1.0:{
                                    'storage_state_binding_method': ['leap_over_forward_only', [2]]
                                }
                            },
                            'boundStartToEnd':{
                                1.0:{
                                    'storage_state_binding_method': ['leap_over_within_solve', [2]] #not sure if exactly the same
                                }
                            }
                        }
                    }
                        #unit__node__profile:
                        #  profile_method:
                        #    upper_limit: 
                        #      profile_method: ['upper_limit', [1, 2, 1]]
                        #    lower_limit: 
                        #      profile_method: ['lower_limit', [1, 2, 1]]
                        #    fixed:
                        #      profile_method: ['fixed', [1, 2, 1]]
                            #no_profile: 
                            #  profile_method: [no_profile, [1, 2, 1]]
                }
            if convertname == 'backbone_to_ines_parameters_to_relationships':
                convertdict = {
                }
            if convertname == 'backbone_to_ines_parameters_to_parameters':
                convertdict = {
                }
            returnlist.append(convertdict)
            if convertdict:
                if save:
                    # save data to a file
                    with open(convertpath,'w') as file:
                        yaml.safe_dump(convertdict, file)
            else:
                print('The file does not exist and neither does the data for ' + convertname)
    return returnlist

def add_base_entities(target_db):
    
    ines_transform.assert_success(target_db.add_alternative_item(name=settings['alternative']), warn=True)
    return target_db

def create_timeline(source_db, target_db):

    timestep_names = range_to_array(settings["timestep_range"])
    sample_names = range_to_array(settings["sample_range"])
    forecast_names = range_to_array(settings["forecast_range"])
    z_names = range_to_array(settings["z_range"])

    #system
    time_value = api.Map(timestep_names,[settings["stepLengthInHours"] for i in timestep_names])
    ines_transform.assert_success(target_db.add_entity_item(entity_class_name='system',entity_byname=('Time',)), warn=True)
    target_db = ines_transform.add_item_to_DB(target_db, "timeline", [settings["alternative"], ("Time",), "system"], time_value, value_type="map")

    #period
    for i in range(0,settings["sample_count"]):
        sample = sample_names[i]
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='period',entity_byname=(sample,)), warn=True)
        target_db = ines_transform.add_item_to_DB(target_db, "start_time", [settings["alternative"], (sample,), "period"], timestep_names[settings["ms_start"][sample]])
        period_duration = (settings["ms_end"][sample]-settings["ms_start"][sample])*settings["stepLengthInHours"]
        target_db = ines_transform.add_item_to_DB(target_db, "duration", [settings["alternative"], (sample,), "period"], period_duration)
        if settings["p_msProbability"][sample]:
            sample_weight = settings["p_msProbability"][sample]
        elif settings["p_msWeight"][sample]:        #check that these actually can be transferred like this
            sample_weight = settings["p_msWeight"][sample]
        elif settings["p_msAnnuityWeight"][sample]:
            sample_weight = settings["p_msAnnuityWeight"][sample]
        else:
            sample_weight = 1
        target_db = ines_transform.add_item_to_DB(target_db, "years_represented", [settings["alternative"], (sample,), "period"], sample_weight)

    #solve_pattern
    ines_transform.assert_success(target_db.add_entity_item(entity_class_name='solve_pattern',entity_byname=('solve',)), warn=True)
    target_db = ines_transform.add_item_to_DB(target_db, "period", [settings["alternative"], ('solve',), "solve_pattern"], api.Array(sample_names[0:settings["sample_count"]]))

    if settings["t_horizon"] and settings["t_jump"]:
        target_db = ines_transform.add_item_to_DB(target_db, "solve_mode", [settings["alternative"], ('solve',), "solve_pattern"], "rolling_window" )
        target_db = ines_transform.add_item_to_DB(target_db, "rolling_jump", [settings["alternative"], ('solve',), "solve_pattern"], settings["t_jump"])
        target_db = ines_transform.add_item_to_DB(target_db, "rolling_horizon", [settings["alternative"], ('solve',), "solve_pattern"], settings["t_horizon"])

    return target_db, timestep_names

def range_to_array(name_range):

    strlength = len(name_range[0])
    letter = name_range[0][0]
    start_str = name_range[0][1:].lstrip('0')
    if start_str == '':
        start = 0
    else:
        start = int(start_str)
    end = int(name_range[1][1:].lstrip('0'))

    name_array = list()
    for i in range(0, end-start):
        num = start + i
        name_array.append(letter+str(num).rjust(strlength-1, '0'))

    return name_array

def create_unit_relationship(source_db, target_db):
    
    params = { # first create relationships
    #'annuity': # not in ines calculate?
    #'availabilityCapacityMargin': Availability of the unit in the capacity margin equation (p.u.). If zero, v_gen is used. Currently used only for output capacity.
    'conversionCoeff': 'conversion_coefficient',
    'fomCosts': 'fixed_cost',
    'inertia': 'inertia_constant',
    #'initialInertia': not in ines hot start
    #'useInitialGeneration': flag for hot start
    'invCosts': 'investment_cost',
    'maxRampDown': 'ramp_limit_down',
    'maxRampUp': 'ramp_limit_up',
    #'startCostCold': 
    #'startCostHot':
    #'startCostWarm'
    #'startFuelConsCold':
    #'startFuelConsHot':
    #'startFuelConsWarm':
    #'unitSizeMVA': not in ines (M volt amp)
    #'upperLimitCapacityRatio': Ratio of the upper limit of the node state and the unit capacity investment ([v_state]/MW)
    'vomCosts': 'other_operational_cost',
    }
    parameter_values = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io')
    for source_entity in source_db.get_entity_items(entity_class_name='grid__node__unit__io'):
        print(source_entity["name"])
        inout = source_entity["entity_byname"][3]
        #alt_ent_class_source = [alt["name"], source_entity["entity_byname"],'grid__node__unit__io']
        if inout == "input":
            target_class = "node__to_unit"
            target_entity_byname = (source_entity["entity_byname"][1],source_entity["entity_byname"][2])
        else:
            target_class = "unit__to_node"
            target_entity_byname = (source_entity["entity_byname"][2],source_entity["entity_byname"][1])

        ines_transform.assert_success(target_db.add_entity_item(entity_class_name=target_class,entity_byname=target_entity_byname), warn=True)
        for source_param_name, target_param_name in params.items():
            for parameter_value in parameter_values:
                if parameter_value["parameter_definition_name"] == source_param_name and parameter_value["entity_name"] == source_entity["name"]:
                    target_db.add_parameter_value_item( entity_class_name=target_class,
                                                            parameter_definition_name=target_param_name,
                                                            entity_byname=target_entity_byname,
                                                            alternative_name=parameter_value["alternative_name"],
                                                            value=parameter_value["value"],
                                                            type=parameter_value["type"])
    
    return target_db

def process_capacities(source_db, target_db):

    alternatives = source_db.get_alternative_items()
    capacities = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='capacity')
    unit_sizes = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='unitSize')
    unit_counts = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='unitCount')
    
    for source_entity in source_db.get_entity_items(entity_class_name='unit'):
        added = False
        for alt in alternatives:
            if not added:
                alt_ent_class_unit = [alt["name"], source_entity["entity_byname"],'unit']
                if not any(source_entity["name"] == unit_count["entity_name"] for unit_count in unit_counts):
                    target_db = ines_transform.add_item_to_DB(target_db, 'units_existing', alt_ent_class_unit, 1, value_type=True)
                    added = True
    
    for source_entity in source_db.get_entity_items(entity_class_name='grid__node__unit__io'):
        if source_entity["entity_byname"][3] == "input":
            target_class_name= "node__to_unit"
            target_entity_byname = (source_entity["entity_byname"][1],source_entity["entity_byname"][2])
        elif source_entity["entity_byname"][3] == "output":
            target_class_name= "unit__to_node"
            target_entity_byname = (source_entity["entity_byname"][2],source_entity["entity_byname"][1])

        unit_capacity = None
        for capacity in capacities:
            if source_entity["name"] == capacity["entity_name"]:
                capacity_value = api.from_database(capacity["value"], capacity["type"])
                unit_count_value=1
                for unit_count in unit_counts:
                    if source_entity["name"] == unit_count["entity_name"]:
                        unit_count_value = api.from_database(unit_count["value"], unit_count["type"])
                unit_capacity =  capacity_value/unit_count_value
                alt_ent_class_target = [capacity["alternative_name"], target_entity_byname, target_class_name]
                target_db = ines_transform.add_item_to_DB(target_db, 'capacity', alt_ent_class_target, unit_capacity, value_type=True)
        if not unit_capacity:
            for unit_size in unit_sizes:
                if source_entity["name"] == unit_size["entity_name"]:
                    target_db.add_parameter_value_item(entity_class_name=alt_ent_class_target[2],
                                                        parameter_definition_name='capacity',
                                                        entity_byname=alt_ent_class_target[1],
                                                        alternative_name=unit_size["alternative_name"],
                                                        value=unit_size["value"],
                                                        type=unit_size["type"])
    return target_db

def create_price_change(source_db, target_db, timesteps):
    alternatives = source_db.get_alternative_items()
    for source_entity in source_db.get_entity_items(entity_class_name='node'):
        for alt in alternatives:
            alt_ent_class_source = [alt["name"], source_entity["entity_byname"],'node']
            alt_ent_class_target = [alt["name"], source_entity["entity_byname"],'node']
            priceChange_map = ines_transform.get_parameter_from_DB(source_db, 'priceChange', alt_ent_class_source)
            if priceChange_map:
                if len(priceChange_map.values) == 1:
                    price = priceChange_map.values[0]
                    target_db = ines_transform.add_item_to_DB(target_db, 'commodity_price', alt_ent_class_target, price, value_type=True)
                else:   #create a timeseries map
                    values = list()
                    price = 0
                    priceChange_dict = priceChange_map.to_dict()
                    for i in timesteps:
                        for step__price in priceChange_dict["data"]:
                            if i == step__price[0]:
                                price = step__price[1]
                                break
                        values.append(price)
                    price_timeseries = api.Map(timesteps, values)
                    target_db = ines_transform.add_item_to_DB(target_db, 'commodity_price', alt_ent_class_target, price_timeseries, value_type="map")
    return target_db

def create_profiles(source_db, target_db):

    alternatives = source_db.get_alternative_items()
    flow_units = source_db.get_entity_items(entity_class_name='flow__unit')
    relationships = source_db.get_entity_items(entity_class_name='grid__node__unit__io')
    fixed_flows = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name = 'fixedFlow')
    for alt in alternatives:
        for flow_node in source_db.get_entity_items(entity_class_name='flow__node'):
            capacity_factor = ines_transform.get_parameter_from_DB(source_db, 'capacityFactor', [alt["name"], flow_node["entity_byname"], 'flow__node'])
            if capacity_factor:
                for flow_unit in flow_units:
                    if flow_node["entity_byname"][0] == flow_unit["entity_byname"][0]:
                        for relationship in relationships:
                            if relationship["entity_byname"][1] == flow_node["entity_byname"][1] and relationship["entity_byname"][2] == flow_unit["entity_byname"][1]:
                                if relationship["entity_byname"][3] == "input":
                                    alt_ent_class_target = [alt["name"], (flow_node["entity_byname"][1], flow_unit["entity_byname"][1]), "node__to_unit"]
                                elif relationship["entity_byname"][3] == "output":
                                    alt_ent_class_target = [alt["name"], (flow_unit["entity_byname"][1], flow_node["entity_byname"][1]), "unit__to_node"]
                                
                                if any(flow_unit["entity_byname"][1] == fixed_flow["entity_name"] for fixed_flow in fixed_flows):
                                    target_db = ines_transform.add_item_to_DB(target_db, 'profile_fix', alt_ent_class_target, capacity_factor, value_type="map")
                                    target_db = ines_transform.add_item_to_DB(target_db, 'profile_method', alt_ent_class_target, "fixed")
                                else:
                                    target_db = ines_transform.add_item_to_DB(target_db, 'profile_limit_upper', alt_ent_class_target, capacity_factor, value_type="map")
                                    target_db = ines_transform.add_item_to_DB(target_db, 'profile_method', alt_ent_class_target, "upper_limit")
    return target_db

def process_links(source_db, target_db):

    #add parameters to entity
    parameters_dict = {    
        # 'annuity' # not in ines calculate?
        'availability': 'availability', 
        'invCost': 'investment_cost',
        #'investMIP': #method for MIP invest
        #'portion_of_transfer_to_reserve' Portion of the import/export that needs to be available as reserve if the interconnection fails
        'variableTransCost': 'operational_cost'
    }

    direct_parameters = {}
    for source_name in parameters_dict.keys():
        direct_parameters[source_name] = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name=source_name)

    transferCaps = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transferCap')
    unit_sizes = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='unitSize')
    transferCapInvLimits = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transferCapInvLimit')
    transfer_loss = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transfer_loss')
    ICrampDown = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='ICrampDown')
    ICrampUp = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='ICrampUp')

    for link in source_db.get_entity_items(entity_class_name='grid__node__node'):
        target_entity_byname = ('link_'+link["entity_byname"][1]+"_"+link["entity_byname"][2],)
        
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='link', entity_byname=target_entity_byname), warn=True)
        for name, target_name in parameters_dict.items():
            for param in direct_parameters[name]:
                if link["name"] == param["entity_name"]:
                    value = api.from_database(param["value"], param["type"])
                    alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
                    target_db = ines_transform.add_item_to_DB(target_db, target_name, alt_ent_class_target, value)

        unit_size_value = None
        for unit_size in unit_sizes:
            if unit_size["entity_name"] == link["name"]:
                unit_size_value = api.from_database(unit_size["value"], unit_size["type"])
                alt_ent_class_target = [unit_size["alternative_name"], target_entity_byname, 'link']
                target_db = ines_transform.add_item_to_DB(target_db, 'capacity', alt_ent_class_target, unit_size_value)
                break
        for transferCap in transferCaps:
            if transferCap["entity_name"] == link["name"]:
                transferCap_value = api.from_database(transferCap["value"], transferCap["type"])
                alt_ent_class_target = [transferCap["alternative_name"], target_entity_byname, 'link']
                if unit_size_value:
                    link_count = transferCap_value/unit_size_value
                    target_db = ines_transform.add_item_to_DB(target_db, 'links_existing', alt_ent_class_target, link_count)
                else:
                    link_count = 1
                    target_db = ines_transform.add_item_to_DB(target_db, 'links_existing', alt_ent_class_target, link_count)
                    target_db = ines_transform.add_item_to_DB(target_db, 'capacity', alt_ent_class_target, transferCap_value)
        for transferCapInvLimit in transferCapInvLimits:
            if transferCapInvLimit["entity_name"] == link["name"]:
                transferCapInvLimit_value = api.from_database(transferCapInvLimit["value"], transferCapInvLimit["type"])
                alt_ent_class_target = [transferCapInvLimit["alternative_name"], target_entity_byname, 'link']
            if unit_size_value:
                links_max_count = (transferCapInvLimit_value + link_count)/unit_size_value
                target_db = ines_transform.add_item_to_DB(target_db, 'links_max_cumulative', alt_ent_class_target, links_max_count)
            else:
                links_max_count = 1
                target_db = ines_transform.add_item_to_DB(target_db, 'links_max_cumulative', alt_ent_class_target, links_max_count)

        for param in transfer_loss:
            if link["name"] == param["entity_name"]:
                value = api.from_database(param["value"], param["type"])
                if isinstance(value,float):
                    value = 1 - value
                if isinstance(value, map):
                    value.values = [1 - val for val in value.values] 
                alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
                target_db = ines_transform.add_item_to_DB(target_db, 'efficiency' , alt_ent_class_target, value)

        for param in ICrampDown:
            if link["name"] == param["entity_name"]:
                value = api.from_database(param["value"], param["type"])
                alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
                target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_down', alt_ent_class_target, 60 * value)

        for param in ICrampUp:
            if link["name"] == param["entity_name"]:
                value = api.from_database(param["value"], param["type"])
                alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
                target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_up', alt_ent_class_target, 60 * value)

    return target_db

def diff_coeff(source_db, target_db):

    for param in source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name = 'diffCoeff'):
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='node__node', 
                                                                entity_byname=(param["entity_byname"][1], param["entity_byname"][2]), warn=True))
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][1], param["entity_byname"][2]), "node_node"]
        value = api.from_database(param["value"], param["type"])
        target_db = ines_transform.add_item_to_DB(target_db, 'diffusion_coefficient', alt_ent_class_target, value)
    
    return target_db

def capacity_margin(source_db, target_db):

    for param in source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'capacityMargin'):
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set', 
                                                                entity_byname=(param["entity_byname"][2],), warn=True))
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__node', 
                                                                entity_byname=(param["entity_byname"][2],param["entity_byname"][2]), warn=True))
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][2],), "set"]
        value = api.from_database(param["value"], param["type"])
        target_db = ines_transform.add_item_to_DB(target_db, 'capacity_margin', alt_ent_class_target, value)

    return target_db

def create_constraints(source_db, target_db):

    for param in source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name = 'boundStateMaxDiff'):
        value = api.from_database(param["value"], param["type"])
        constraint_map = {
            "less_than": ['state_leq',param["entity_byname"][1], param["entity_byname"][2]],
            "greater_than": ['state_geq', param["entity_byname"][2], param["entity_byname"][1]]
        }
        for name, value in constraint_map:
            target_byname = ("_".join(value[0],value[1],value[2]),)
            ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', 
                                                                    entity_byname=target_byname, warn=True))
            alt_ent_class_constraint = [param["alternative_name"], target_byname, "constraint"]
            target_db = ines_transform.add_item_to_DB(target_db, 'sense', alt_ent_class_constraint, name)
            target_db = ines_transform.add_item_to_DB(target_db, 'constant', alt_ent_class_constraint, value)
            
            node_param = api.Map([target_byname], [1])
            alt_ent_class_node = [param["alternative_name"],(param["entity_byname"][1],), "node"]
            target_db = ines_transform.add_item_to_DB(target_db, 'constraint_storage_state_coefficient', alt_ent_class_node, node_param)
            node_param = api.Map([target_byname], [-1])
            alt_ent_class_node = [param["alternative_name"],(param["entity_byname"][2],), "node"]
            target_db = ines_transform.add_item_to_DB(target_db, 'constraint_storage_state_coefficient', alt_ent_class_node, node_param)
    
    #for param in source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_value_name = 'upperLimitCapacityRatio'):

    return target_db

def add_node_types(source_db, target_db):

    grid_node_boudaries =  source_db.get_entity_items(entity_class_name='grid__node__boundary')
    prices = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'price')
    price_changes = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'priceChange')
    node_balances = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'nodeBalance')
    use_constants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useConstant')
    use_timeseries = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useTimeSeries')
    

    for node in source_db.get_entity_items(entity_class_name='node'):
        node_type = None
        found = False
        for price in prices:
            if price["entity_name"]== node["name"]:
                node_type = "commodity"
                alt = price["alternative_name"]
                found = True
                break
        if not found:
            for price_change in price_changes:
                if price_change["entity_name"]== node["name"]:
                    node_type = "commodity"
                    alt = price_change["alternative_name"]
                    found = True
                    break
        if not found:
            for grid_node_boudary in grid_node_boudaries:
                if ((grid_node_boudary["entity_byname"][2] == 'upwardLimit' or 
                        grid_node_boudary["entity_byname"][2] == 'downwardLimit') and
                        grid_node_boudary["entity_byname"][1] == node["name"]):
                    for use_constant in use_constants:
                        if use_constant["entity_name"] == grid_node_boudary["name"] and api.from_database(use_constant["value"],use_constant["type"]) == 1.0:
                            node_type = "storage" 
                            alt =  use_constant["alternative_name"]
                            found = True
                            break
                    if not found:
                        for use_timeserie in use_timeseries:
                            if use_timeserie["entity_name"] == grid_node_boudary["name"] and api.from_database(use_timeserie["value"],use_timeserie["type"]) == 1.0:
                                node_type = "storage" 
                                alt = use_timeserie["alternative_name"]
                                found = True
                                break  
        if not found:
            for node_balance in node_balances:
                if node_balance["entity_byname"][1]== node["name"] and api.from_database(node_balance["value"],node_balance["type"]) == 1.0:
                    node_type = "balance"
                    alt = node_balance["alternative_name"]
                    break

        target_db = ines_transform.add_item_to_DB(target_db, 'node_type', [alt, node["entity_byname"], "node"], node_type)
    
    return target_db

if __name__ == "__main__":
    developer_mode = False
    if developer_mode:
        # save entities to yaml file
        save_folder = os.path.dirname(__file__)
        conversion_configuration(conversions = [save_folder+'/backbone_to_ines_entities.yaml', save_folder+'/backbone_to_ines_parameters.yaml', save_folder+'/backbone_to_ines_parameter_methods.yaml',
                                             save_folder+'/backbone_to_ines_parameters_to_relationships.yaml'], overwrite=True)
    else:
        # assume the file to be used inside of Spine Toolbox
        #url_db_in = sys.argv[1]
        #url_db_out = sys.argv[2]
        #settings_path = 'backbone_to_ines_settings.yaml'

        url_db_in = 'sqlite:///C:/Users/aetart/Documents/ines-backbone/BB_data_test_debug.sqlite'
        url_db_out = 'sqlite:///C:/Users/aetart/Documents/ines-backbone/ines-spec.sqlite'
        settings_path = 'C:/Users/aetart/Documents/ines-backbone/backbone-ines/backbone_to_ines_settings.yaml'

        # open yaml files
        entities_to_copy,parameter_transforms,parameter_methods, parameters_to_relationships, parameters_to_parameters = conversion_configuration()
        with open(settings_path,'r') as file:
            settings = yaml.safe_load(file) 

        main()