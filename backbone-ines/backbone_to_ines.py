import os
import sys
import yaml
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from ines_tools import ines_transform
from ines_tools import ines_initialize
import spinedb_api as api
from sqlalchemy.exc import DBAPIError
from spinedb_api.exception import NothingToCommit

def main():
    # transform spine db with backbone data (source db) into a spine db that already has the ines structure (target_db)
    with api.DatabaseMapping(url_db_in) as source_db:
        with api.DatabaseMapping(url_db_out) as target_db:
            # completely empty database
            #purge.purge(target_db, purge_settings=None)
            # add ines structure
            # empty database except for ines structure
            target_db = ines_initialize.purge_db_from_data(target_db)
            source_db = ines_initialize.fetch_data(source_db)
            target_db = ines_initialize.copy_alternatives_scenarios(source_db, target_db)

            #Add entities that are always there
            target_db = add_base_entities(target_db)

            target_db, t_val__timestamp = create_timeline(source_db, target_db)

            # copy entities from yaml files
            print("copy entities")
            target_db = ines_transform.copy_entities(source_db, target_db, entities_to_copy)

            #target_db = inflow_timeseries_from_csv(target_db, "./influx_timeseries", t_val__timestamp)
            #create in-and-out relationships
            print("add link capacities")
            target_db = process_links(source_db, target_db, t_val__timestamp)
            print("create unit relationships")
            target_db = create_unit_relationship(source_db, target_db, t_val__timestamp)
            print("create grid groups")
            target_db = create_sets_from_grids(source_db, target_db)
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
            print("unit parameters")
            target_db = create_unit_parameters(source_db, target_db, t_val__timestamp)
            print("add capacities")
            target_db = process_capacities(source_db, target_db, t_val__timestamp)
            target_db = create_node_capacities(source_db, target_db, t_val__timestamp)
            print("create timeseries from price change")
            target_db = create_price_change(source_db, target_db, t_val__timestamp)
            print("create profiles")
            target_db = create_profiles(source_db, target_db, t_val__timestamp)
            print("create emissions")
            target_db = create_emissions(source_db, target_db, t_val__timestamp)
            print("create reserves")
            target_db = create_reserves(source_db, target_db, t_val__timestamp)
            print("add capacity margin")
            target_db = capacity_margin(source_db, target_db, t_val__timestamp)
            print("create user constraints")
            target_db = create_unit_node_constraints(source_db, target_db, t_val__timestamp)
            target_db = create_bound_state_constraints(source_db, target_db, t_val__timestamp)
            print("handle boundaries")
            target_db = handle_boundaries(source_db, target_db, t_val__timestamp)
            print("diffusion coefficient")
            target_db = diff_coeff(source_db, target_db, t_val__timestamp)
            print("group constraints")
            target_db = create_group_constraints(source_db, target_db)
            print("add simple timeseries")
            target_db = create_simple_timeseries(source_db, target_db, t_val__timestamp)
            print("add node types")
            target_db = add_node_types(source_db, target_db)
            print("add entity alternative items")
            target_db = add_entity_alternative_items(target_db)
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

            if len(cannot_convert) > 0:
                print(cannot_convert)
            print("transformation complete")

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
                    'group': ['set'],
                    'restype': ['reserve'],
                    'unit': ['unit'],
                    'node': ['node'],
                    'grid__node__group': {'set__node': [[3], [2]]},
                    'unit__group': {'set__unit': [[2], [1]]}
                }
            if convertname == 'backbone_to_ines_parameters':
                convertdict = {
                    'constraint': {
                    },
                    'group':{
                        'set':{
                            #ROCOF:  #inertia stuff
                            #'constrainedCapTotalMax': 'invest_max_total', #check units
                            #'defaultFrequency':  inertia stuff
                            #'dynamicInertia' #flag, not in ines
                            'energyMax': 'flow_max_cumulative',
                            'energyMin':'flow_min_cumulative', 
                            #'energyShareMax'
                            #'energyShareMin'
                            #'staticInertia' flag, part of the ROCOF stuff
                        }
                    },
                    'model':{ 
                        #'system':{
                            #'discountFactor': not 'discount_rate'! calculate discount rate? period wise discount rate??
                            #'t_invest': #assumed to be the period start in ines
                        #}
                    },
                    'node':{
                    },
                    'restype':{
                        #'restype':{
                            #'restypeReleasedForRealization'
                            #'restype_inertia'
                        #}
                    },
                    'unit':{
                        'unit':{
                            #'becomeAvailable':  not ines or calculate to the availability timeseries
                            #'becomeUnavailable': not ines or calculate to the availability timeseries
                            #'boundSamples': not in ines
                            #'initialOnlineStatus': #hot start not in ines
                            'maxUnitCount': 'units_max_cumulative',
                            'minOperationHours': ['min_uptime', 60],
                            'minShutdownHours': ['min_downtime', 60],
                            'minUnitCount': 'units_min_cumulative',
                            #'rampSpeedFromMinLoad': not in ines
                            #'rampSpeedToMinLoad' not in ines
                            'unitCount': 'units_existing',
                            #'useInitialOnlineStatus' #ines hot start flag
                        }
                    },
                    'grid__node':{
                        'node':{
                            'selfDischargeLoss': ['storage_loss_from_stored_energy', 1.0, [[2]]]
                        }
                    },
                    'group__emission':{
                    },
                    'group__restype':{
                        #update_offset
                    },
                    'node__emission':{
                    },
                    'unit__group':{
                        #'constrainedCapMultiplier': #for constrainedCapTotalMax, not in ines 
                    },
                    'unit__node':{
                    },
                    'grid__node__boundary':{
                        #slackCost
                    },
                    'grid__node__group':{
                    },
                    'grid__restype__up_down':{
                    },
                    'grid__node__unit__boundary':{
                    },
                    'grid__node__unit__emission':{
                        #'invEmissionFactor':   investment emissions possibility to divide over years
                    },
                    'grid__node__unit__restype':
                        {
                            #offlineReserveCapability
                        }
                    ,'grid__node__unit__io':{
                        'unit':{
                            'shutdownCost': ['shutdown_cost',1.0,[[3]]]
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
                                    'storage_state_fix_method': ['fix_start', [2]]
                                }
                            },
                            'boundStartAndEnd':{
                                1.0:{
                                    'storage_state_fix_method': ['fix_start_and_horizon_end', [2]]
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

    #timestamp map
    t_val__timestamp = dict()
    timestamps = list()
    if settings["t_start"] and settings["t_end"]:
        for i, step in enumerate(timestep_names):
            if i >= settings["t_start"] and i < settings["t_end"]:
                timestamp = datetime.fromisoformat(settings["first_timestamp"]) +  timedelta(hours=i*settings["stepLengthInHours"])
                timestamps.append(timestamp)
                t_val__timestamp[step] = timestamp
    elif settings["tsYear"] and settings["modelledDays"]:
        start_timestep = datetime(settings["tsYear"],1,1)
        start_time = (start_timestep - datetime.fromisoformat(settings["first_timestamp"])).total_seconds()/3600.0
        end_timestep = start_timestep + timedelta(hours=settings["modelledDays"]*settings["stepLengthInHours"])
        end_time = (end_timestep - datetime.fromisoformat(settings["first_timestamp"])).total_seconds()/3600.0
        for i, step in enumerate(timestep_names):
            if i>= start_time and i < end_time:
                timestamp = start_timestep +  timedelta(hours=(i-start_time)*settings["stepLengthInHours"])
                timestamps.append(timestamp)
                t_val__timestamp[step] = timestamp
    else:
        print("Define t_start and t_end in the settings file. Alternatively, if using Northern European Model, set tsYear and modelledDays")
        exit(-1)
    #system
    time_series = api.TimeSeriesVariableResolution(timestamps, [settings["stepLengthInHours"] for i in timestamps], ignore_year = False, repeat=False, index_name="time step")
    ines_transform.assert_success(target_db.add_entity_item(entity_class_name='system',entity_byname=('Time',)), warn=True)
    target_db = ines_transform.add_item_to_DB(target_db, "timeline", [settings["alternative"], ("Time",), "system"], time_series, value_type="map")

    #period
    start_times = []
    period_durations = []
    for i in range(0,settings["sample_count"]):
        sample = sample_names[i]
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='period',entity_byname=(sample,)), warn=True)
        period_duration = api.Duration(str((settings["ms_end"][sample]-settings["ms_start"][sample])*settings["stepLengthInHours"])+"h")
        start_time = api.DateTime(timestamps[settings["ms_start"][sample]-1])
        start_times.append(start_time)
        period_durations.append(period_duration)
        if settings["p_msProbability"][sample]:
            sample_weight = settings["p_msProbability"][sample]
        elif settings["p_msWeight"][sample]:        #check that these actually can be transferred like this
            sample_weight = settings["p_msWeight"][sample]
        elif settings["p_msAnnuityWeight"][sample]:
            sample_weight = settings["p_msAnnuityWeight"][sample]
        else:
            sample_weight = 1
        target_db = ines_transform.add_item_to_DB(target_db, "years_represented", [settings["alternative"], (sample,), "period"], sample_weight)
        target_db = ines_transform.add_item_to_DB(target_db, "start_time", [settings["alternative"], (sample,), "period"], start_time)
    #solve_pattern
    ines_transform.assert_success(target_db.add_entity_item(entity_class_name='solve_pattern',entity_byname=('solve',)), warn=True)
    target_db = ines_transform.add_item_to_DB(target_db, "period", [settings["alternative"], ('solve',), "solve_pattern"], api.Array(sample_names[0:settings["sample_count"]]))
    target_db = ines_transform.add_item_to_DB(target_db, "start_time", [settings["alternative"], ('solve',), "solve_pattern"], api.Array(start_times))
    target_db = ines_transform.add_item_to_DB(target_db, "duration", [settings["alternative"], ('solve',), "solve_pattern"], api.Array(period_durations))

    if settings["t_horizon"] and settings["t_jump"]:
        target_db = ines_transform.add_item_to_DB(target_db, "solve_mode", [settings["alternative"], ('solve',), "solve_pattern"], "rolling_window")
        jump = api.Duration(str(settings["t_jump"] * settings["stepLengthInHours"])+"h")
        horizon = api.Duration(str(settings["t_horizon"] * settings["stepLengthInHours"])+"h")
        target_db = ines_transform.add_item_to_DB(target_db, "rolling_jump", [settings["alternative"], ('solve',), "solve_pattern"], jump)
        target_db = ines_transform.add_item_to_DB(target_db, "rolling_horizon", [settings["alternative"], ('solve',), "solve_pattern"], horizon)
    if settings["stepLengthInHours"]:
        time_resolution = api.Duration(str(settings["stepLengthInHours"])+"h")
        target_db = ines_transform.add_item_to_DB(target_db, "time_resolution", [settings["alternative"], ('solve',), "solve_pattern"], time_resolution)
    #stochastic information
    ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set',entity_byname=('stochastics',)), warn=True)
    target_db = ines_transform.add_item_to_DB(target_db, "stochastic_forecast_weights", [settings["alternative"], ('stochastics',), "set"],
                                              api.Map(list(settings["p_mfProbability"].keys()),list(settings["p_mfProbability"].values())))
    target_db = ines_transform.add_item_to_DB(target_db, "stochastic_scope", [settings["alternative"], ('solve',), "solve_pattern"], "whole_model")

    #forecast interpolation
    if settings["t_improveForecastNew"]:
        inter_array = api.Array([i/float(settings["t_improveForecastNew"]) for i in range(0,settings["t_improveForecastNew"])])
        target_db = ines_transform.add_item_to_DB(target_db, "stochastic_forecast_interpolation_factors", [settings["alternative"], ('stochastics',), "set"], inter_array)
        target_db = ines_transform.add_item_to_DB(target_db, "stochastic_method", [settings["alternative"], ('stochastics',), "set"], "interpolate_time_series_forecasts")
        
    
    return target_db, t_val__timestamp

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

def create_unit_relationship(source_db, target_db, t_val__timestamp):
    
    params = { # first create relationships
    #'availabilityCapacityMargin': Availability of the unit in the capacity margin equation (p.u.). If zero, v_gen is used. Currently used only for output capacity.
    'conversionCoeff': 'conversion_coefficient',
    'fomCosts': 'fixed_cost',
    'inertia': 'inertia_constant',
    #'initialInertia': not in ines hot start
    #'initialGeneration': not in ines hot start
    #'useInitialGeneration': flag for hot start
    'invCosts': 'investment_cost',
    'maxRampDown': 'ramp_limit_down',
    'maxRampUp': 'ramp_limit_up',
    #'unitSizeMVA': not in ines (M volt amp)
    }

    grid__node__unit__groups = source_db.get_entity_items(entity_class_name='grid__node__unit_group')
    annuity = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='annuity')
    parameter_values = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io')
    units = source_db.get_entity_items(entity_class_name='unit')
    gnuios = source_db.get_entity_items(entity_class_name='grid__node__unit__io')
    vomCosts = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='vomCosts') 
    priceChanges = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name='priceChange')
    emissionContents = source_db.get_parameter_value_items(entity_class_name='node__emission', parameter_definition_name='emission_content')
    fixedFuelFractions = source_db.get_parameter_value_items(entity_class_name='unit__node', parameter_definition_name='fixedFuelFractions')

    startCostCold = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startCostCold')
    startCostHot = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startCostHot') 
    startCostWarm = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startCostWarm') 
    startFuelConsCold = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startFuelConsCold')
    startFuelConsHot = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startFuelConsHot') 
    startFuelConsWarm = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='startFuelConsWarm') 

    for source_entity in gnuios:
        print(source_entity["name"])
        inout = source_entity["entity_byname"][3]
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
        for param in vomCosts:
            if param["entity_byname"] == source_entity["entity_byname"]:
                value = api.from_database(param["value"], param["type"])
                target_db = pass_timeseries(target_db,'other_operational_cost', 'other_operational_cost_forecasts', value, 
                                            [param["alternative_name"],target_entity_byname,target_class], t_val__timestamp)
                
        for param in annuity:
            if param["entity_byname"] == source_entity["entity_byname"]:
                value = api.from_database(param["value"], param["type"])
                target_db = calculate_investment_cost(source_db, target_db, [param["alternative_name"], param["entity_byname"], target_class], value, storage = False)
                lifetime = settings["default_lifetime"]
                r = settings["default_interest_rate"]
                target_db = ines_transform.add_item_to_DB(target_db, 'lifetime', [param["alternative_name"], (param["entity_byname"][2],), target_class], lifetime)
                target_db = ines_transform.add_item_to_DB(target_db, 'interest_rate', [param["alternative_name"], (param["entity_byname"][2],), target_class], r)
    
        for entity in grid__node__unit__groups:
            group_entity_byname = (entity["entity_byname"][3],target_entity_byname[0], target_entity_byname[1])
            ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__unit_flow', 
                                                entity_byname=group_entity_byname), warn=True)

    #calculate start-up costs and emissions
    for unit in units:
        emission_for_fuel = dict() 
        price_for_fuel = dict()
        #get inputfuels
        for source_entity in gnuios:
            if unit["entity_byname"][0] == source_entity["entity_byname"][2] and source_entity["entity_byname"][3] == "input":
                fixed_fuel_fraction = 1
                for fixedFuelFraction in fixedFuelFractions:
                    if source_entity["entity_byname"][2] == fixedFuelFraction["entity_byname"][0] and source_entity["entity_byname"][1] == fixedFuelFraction["entity_byname"][1]:
                        fixed_fuel_fraction = api.from_database(fixedFuelFraction["value"], fixedFuelFraction["type"])
                for emissionContent in emissionContents:
                    if source_entity["entity_byname"][1] == emissionContent["entity_byname"][0] and emissionContent["entity_byname"][1] == "CO2":
                        value = api.from_database(emissionContent["value"], emissionContent["type"])
                        emission_for_fuel[emissionContent["entity_byname"][0]] = value * fixed_fuel_fraction

                for priceChange in priceChanges:
                    if source_entity["entity_byname"][1] == priceChange["entity_byname"][0]:
                        #only the first value is taken
                        #will not work if a true time series
                        #would need a change in ines-spec
                        value = api.from_database(priceChange["value"], priceChange["type"])
                        price_for_fuel[priceChange["entity_byname"][0]] = float(value.values[0]) * fixed_fuel_fraction
        
        emission_for_fuels = sum(emission_for_fuel.values()) 
        price_for_fuels = sum(price_for_fuel.values())

        emissions = dict()
        prices = dict()
        for source_entity in gnuios:
            if unit["entity_byname"][0] == source_entity["entity_byname"][2]:
                start_up_tiers = dict()
                start_up_fuel_tiers = dict()
                for start_hot in startCostHot:
                    if start_hot["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_hot["value"], start_hot["type"])
                        start_up_tiers[0] =  value
                        alt = start_hot["alternative_name"]
                for start_warm in startCostWarm:
                    if start_warm["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_warm["value"], start_warm["type"])
                        start_up_tiers[1] =  value
                for start_cold in startCostCold:
                    if start_cold["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_cold["value"], start_cold["type"])
                        start_up_tiers[2] =  value
                        alt = start_cold["alternative_name"]
                for start_fuel_hot in startFuelConsHot:
                    if start_fuel_hot["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_fuel_hot["value"], start_fuel_hot["type"])
                        start_up_fuel_tiers[0] =  value
                        alt = start_cold["alternative_name"]
                for start_fuel_warm in startFuelConsWarm:
                    if start_fuel_warm["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_fuel_warm["value"], start_fuel_warm["type"])
                        start_up_fuel_tiers[1] =  value
                for start_fuel_cold in startFuelConsCold:
                    if start_fuel_cold["entity_byname"] == source_entity["entity_byname"]:
                        value = api.from_database(start_fuel_cold["value"], start_fuel_cold["type"])
                        start_up_fuel_tiers[2] =  value
                        alt = start_cold["alternative_name"]
                if len(start_up_fuel_tiers.keys()) > 0 or len(start_up_tiers.keys()) > 0:
                    for key, value in start_up_fuel_tiers.items():
                        if key not in emissions.keys():
                            emissions[key] = value * emission_for_fuels
                        else:
                            emissions[key] = emissions[key] + value * emission_for_fuels
                        if key not in prices.keys():
                            prices[key] = value * price_for_fuels
                        else:
                            prices[key] = prices[key] + value * price_for_fuels
                    for key, value in start_up_tiers.items():
                        if key not in prices.keys():
                            prices[key] = value
                        else:
                            prices[key] = prices[key] + value
        if len(list(prices.keys())) > 1:
            out = api.Map([str(x) for x in prices.keys()],list(prices.values()))
            target_db = ines_transform.add_item_to_DB(target_db, 'startup_cost_tiers', [alt, unit["entity_byname"] ,"unit"], out, value_type=True)
        elif len(list(prices.keys())) == 1:
            target_db = ines_transform.add_item_to_DB(target_db, 'startup_cost', [alt, unit["entity_byname"] ,"unit"], prices.popitem()[1], value_type=True)
        if len(list(emissions.keys())) > 1:
            out = api.Map([str(x) for x in emissions.keys()],list(emissions.values()))
            target_db = ines_transform.add_item_to_DB(target_db, 'startup_co2_emission_tiers', [alt, unit["entity_byname"] ,"unit"], out, value_type=True)
        elif len(list(emissions.keys())) == 1:
            target_db = ines_transform.add_item_to_DB(target_db, 'startup_co2_emission', [alt, unit["entity_byname"] ,"unit"], emissions.popitem()[1], value_type=True)
                        
    return target_db

def process_capacities(source_db, target_db, t_val__timestamp):

    alternatives = source_db.get_alternative_items()
    capacities = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name='capacity')
    unit_sizes = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='unitSize')
    unit_counts = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='unitCount')
    
    for source_entity in source_db.get_entity_items(entity_class_name='unit'):
        # if unit count not given, assume 1, unit count is moved otherwise by the helper function
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

        #add the capacity and unit count
        unit_capacity = None
        for capacity in capacities:
            if source_entity["name"] == capacity["entity_name"]:
                capacity_value = api.from_database(capacity["value"], capacity["type"])
                unit_count_value = 1
                for unit_count in unit_counts:
                    if source_entity["name"] == unit_count["entity_name"]:
                        unit_count_value = api.from_database(unit_count["value"], unit_count["type"])
                unit_capacity =  capacity_value/unit_count_value
                alt_ent_class_target = [capacity["alternative_name"], target_entity_byname, target_class_name]
                target_db = ines_transform.add_item_to_DB(target_db, 'capacity', alt_ent_class_target, unit_capacity, value_type=True)
        #if capacity not given, check if unit size is given and use that
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

def create_price_change(source_db, target_db, t_val__timestamp):
    for param in source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'priceChange'):
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"],'node']
        target_db = single_price_change(target_db, t_val__timestamp, value, alt_ent_class_target, 'commodity_price', 'commodity_price_forecasts')
     
    return target_db

def create_profiles(source_db, target_db, t_val__timestamp):

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
                                    target_db = pass_timeseries(target_db,'profile_limit_upper', 'profile_limit_upper_forecasts', capacity_factor, alt_ent_class_target, t_val__timestamp)
                                    target_db = ines_transform.add_item_to_DB(target_db, 'profile_method', alt_ent_class_target, "upper_limit")
    return target_db

def create_unit_parameters(source_db, target_db, t_val__timestamp):

    units = source_db.get_entity_items(entity_class_name='unit')
    investMIP = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='investMIP')
    efficiency  = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='efficiency')
    efficiency_ts  = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name='efficiency_ts')
    unit_availabilities = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name = 'availability')
    startColdAfterXhourss = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name = 'startColdAfterXhours')
    startWarmAfterXHourss = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name = 'startWarmAfterXHours')
    effLevel__effSelector__unit =  source_db.get_entity_items(entity_class_name='effLevel__effSelector__unit')

    for unit in units:
        #first check if tiered unit
        cool_down_tiers = dict()
        tiered = False
        for param in startWarmAfterXHourss:
            if unit["entity_byname"] == param["entity_byname"]:
                cool_down_tiers["1"] = api.from_database(param["value"], param["type"]) * 60
                alt = param["alternative_name"]
        for param in startColdAfterXhourss:
            if unit["entity_byname"] == param["entity_byname"]:
                cool_down_tiers["2"] = api.from_database(param["value"], param["type"]) * 60
                alt = param["alternative_name"]
        if len(list(cool_down_tiers.keys())) > 0:
            val = api.Map(list(cool_down_tiers.keys()),list(cool_down_tiers.values()))
            target_db = ines_transform.add_item_to_DB(target_db, 'cooling_time_to_tiers', [alt, unit["entity_byname"], "unit"], val)
            tiered = True
        
        #Check startup method
        # Here we are using only the effLevel level1 and ignoring the levels used with longer step sizes
        # when rampToMinLoad and rampFromMinLoad are added, these units would use the "trajectory" method
        for entity in effLevel__effSelector__unit:
            if unit["entity_byname"][0] == entity["entity_byname"][2]:
                if entity["entity_byname"][0] == "level1":
                    if entity["entity_byname"][1] == "directOff":
                        method = "no_startup"
                    elif entity["entity_byname"][1] == "directOnLP":
                        if tiered:
                            method = "linear_with_tiers"
                        else:
                            method = "linear"
                    elif entity["entity_byname"][1] == "directOnMIP":
                        if tiered:
                            method = "integer_with_tiers"
                        else:
                            method = "integer"
                    target_db = ines_transform.add_item_to_DB(target_db, 'startup_method', [settings["alternative"], (entity["entity_byname"][2],) ,"unit"], method)
                    break

        for param in efficiency:
            if unit["entity_byname"] == param["entity_byname"]:
                value = api.from_database(param["value"], param["type"])
                effs = []
                ops = []
                for i in range(0,len(value.values)):
                    if value.indexes[i][0:3] == "eff":
                        effs.append(value.values[i])
                    if value.indexes[i][0:2] == "op": 
                        ops.append(value.values[i]) 
                    if value.indexes[i] == "section":
                        ops.insert(0,0)
                        effs.insert(0,value.values[i])

                    if value.indexes[i][0:3] == "hr":
                        effs.append(1/(value.values[i]/3.6))
                    if value.indexes[i][0:2] == "hrop": 
                        ops.append(value.values[i]) 
                    if value.indexes[i] == "hrsection":
                        ops.insert(0,0)
                        effs.insert(0,1/(value.values[i]/3.6))

                    if value.indexes[i] == "opFirstCross":
                        cannot_convert.append("Efficiency parameter opFirstCross cannot be converted in " + param["entity_byname"] +"\n")
                #when single eff value back bone has two eff points and one op
                if method == "no_startup" or len(ops) == 1:
                    outval = max([float(eff) for eff in effs])
                    conversion_method = "constant_efficiency"

                elif len(effs) > 0 and len(ops) > 0:
                    outval = api.Map(ops,effs)
                    if len(ops) == 2:
                        conversion_method = "partial_load_efficiency"
                    else:
                        conversion_method = "piecewise_linear"

                #out
                alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "unit"]
                target_db = ines_transform.add_item_to_DB(target_db, 'efficiency', alt_ent_class_target, outval)
                target_db = ines_transform.add_item_to_DB(target_db, 'conversion_method', alt_ent_class_target, conversion_method)


    #eff00, f, t, val 
    #only eff00 used, constant for each time step
    for param in efficiency_ts:
        value = api.from_database(param["value"], param["type"])
        if len(value.indexes) > 1:
            cannot_convert.append("INES does not support piecewise linear timeseries for efficiency (4d map)")
        else:
            target_db = pass_timeseries(target_db, 'efficiency', 'efficiency_forecasts', value.values, [param["alternative_name"], (param["entity_byname"][0],), "unit"], t_val__timestamp)
   
    for param in unit_availabilities:
        value =  api.from_database(param["value"], param["type"])
        target_db = pass_timeseries(target_db, 'availability', 'availability_forecasts', value, [param["alternative_name"], (param["entity_byname"][0],), "unit"], t_val__timestamp)

    for param in investMIP:
        value = api.from_database(param["value"], param["type"])
        if value != 0.0:
            alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "unit"]
            target_db = ines_transform.add_item_to_DB(target_db, 'investment_uses_integer', alt_ent_class_target, True)

    return target_db

def process_links(source_db, target_db, t_val__timestamp):

    #add parameters to entity
    parameters_dict = {    
        'invCost': 'investment_cost',
    }

    direct_parameters = {}
    for source_name in parameters_dict.keys():
        direct_parameters[source_name] = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name=source_name)

    grid__node__node__groups = source_db.get_entity_items(entity_class_name='grid__node__node__group')

    annuity = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='annuity')
    availabilities = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='availability')
    transferCaps = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transferCap')
    unit_sizes = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='unitSize')
    transferCapInvLimits = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transferCapInvLimit')
    transfer_loss = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='transferLoss')
    ICrampDown = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='ICrampDown')
    ICrampUp = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='ICrampUp')
    investMIP = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='investMIP')
    variableTransCost = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name='variableTransCost')

    for link in source_db.get_entity_items(entity_class_name='grid__node__node'):
        target_entity_byname = ('link_'+link["entity_byname"][1]+"_"+link["entity_byname"][2],)
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='link', entity_byname=target_entity_byname), warn=True)
        rel_target_entity_byname = (link["entity_byname"][1], target_entity_byname[0], link["entity_byname"][2])
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='node__link__node', entity_byname=rel_target_entity_byname), warn=True)
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

    for param in variableTransCost:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value = api.from_database(param["value"], param["type"])
        target_db = pass_timeseries(target_db,'operational_cost', 'operational_cost_forecasts', value, 
                                    [param["alternative_name"], target_entity_byname, "link"], t_val__timestamp)

    for param in transfer_loss:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value = api.from_database(param["value"], param["type"])
        if isinstance(value,float):
            value = 1 - value
        if isinstance(value, map):
            value.values = [1 - val for val in value.values] 
        alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
        target_db = pass_timeseries(target_db, 'efficiency', 'efficiency_forecasts', value, alt_ent_class_target, t_val__timestamp)

    for param in ICrampDown:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
        target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_down', alt_ent_class_target, 60 * value)

    for param in ICrampUp:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
        target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_up', alt_ent_class_target, 60 * value)

    for param in investMIP:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value = api.from_database(param["value"], param["type"])
        if value != 0.0:
            alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
            target_db = ines_transform.add_item_to_DB(target_db, 'investment_uses_integer', alt_ent_class_target, True)
    for param in availabilities:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value =  api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
        target_db = pass_timeseries(target_db, 'availability', 'availability_forecasts', value, alt_ent_class_target, t_val__timestamp)

    for param in annuity:
        target_entity_byname = ('link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2],)
        value =  api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], target_entity_byname, "link"]
        lifetime = settings["default_lifetime"]
        r = settings["default_interest_rate"]
        target_db = ines_transform.add_item_to_DB(target_db, 'lifetime', alt_ent_class_target, lifetime)
        target_db = ines_transform.add_item_to_DB(target_db, 'interest_rate', alt_ent_class_target, r)
        target_db = calculate_investment_cost(source_db, target_db, value, alt_ent_class_target)
    
    for entity in grid__node__node__groups:
        target_entity_byname = (entity["entity_byname"][3],'link_'+entity["entity_byname"][1]+"_"+entity["entity_byname"][2])
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__link', 
                                            entity_byname=target_entity_byname), warn=True)

    return target_db

def diff_coeff(source_db, target_db, t_val__timestamp):

    for param in source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name = 'diffCoeff'):
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='node__node', 
                                                                entity_byname=(param["entity_byname"][1], param["entity_byname"][2])), warn=True)
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][1], param["entity_byname"][2]), "node_node"]
        value = api.from_database(param["value"], param["type"])
        target_db = ines_transform.add_item_to_DB(target_db, 'diffusion_coefficient', alt_ent_class_target, value)
    
    return target_db

def capacity_margin(source_db, target_db, t_val__timestamp):

    for param in source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'capacityMargin'):
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set', 
                                                                entity_byname=(param["entity_byname"][2],)), warn=True)
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__node', 
                                                                entity_byname=(param["entity_byname"][2],param["entity_byname"][2])), warn=True)
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][2],), "set"]
        value = api.from_database(param["value"], param["type"])
        target_db = ines_transform.add_item_to_DB(target_db, 'capacity_margin', alt_ent_class_target, value)
        target_db = ines_transform.add_item_to_DB(target_db, 'capacity_margin_method', alt_ent_class_target, 'dispatch')

    return target_db

def create_emissions(source_db, target_db, t_val__timestamp):

    emissionCaps =  source_db.get_parameter_value_items(entity_class_name='group__emission', parameter_definition_name = 'emissionCap')
    emissionTaxs =  source_db.get_parameter_value_items(entity_class_name='group__emission', parameter_definition_name = 'emissionTax')
    emissionPriceChanges = source_db.get_parameter_value_items(entity_class_name='group__emission', parameter_definition_name = 'emissionPriceChange')
    emission_contents = source_db.get_parameter_value_items(entity_class_name='node__emission', parameter_definition_name = 'emission_content')
    fomEmissions = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__emission', parameter_definition_name = 'emissionCap')
    #invEmissionFactor = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__emission', parameter_definition_name = 'emissionCap')
    invEmissions = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__emission', parameter_definition_name = 'invEmission')
    vomEmissions = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__emission', parameter_definition_name = 'vomEmission')
    grid__node__unit__io = source_db.get_entity_items(entity_class_name='grid__node__unit__io')


    for emission_content in emission_contents:
        value = api.from_database(emission_content["value"], emission_content["type"])
        alt_ent_class_target = [emission_content["alternative_name"],(emission_content["entity_byname"][0],), "node"]
        if emission_content["entity_byname"][1] == "CO2" or emission_content["entity_byname"][1] == "co2":
            target_db = ines_transform.add_item_to_DB(target_db, "co2_content", alt_ent_class_target, value)

    for emissionCap in emissionCaps:
        value = api.from_database(emissionCap["value"], emissionCap["type"])
        alt_ent_class_target = [emissionCap["alternative_name"],(emissionCap["entity_byname"][0],), "group"]
        if emissionCap["entity_byname"][1] == "CO2" or emissionCap["entity_byname"][1] == "co2":
            target_db = ines_transform.add_item_to_DB(target_db, "co2_max_cumulative", alt_ent_class_target, value)
        if emissionCap["entity_byname"][1] == "SO2" or emissionCap["entity_byname"][1] == "so2":
            target_db = ines_transform.add_item_to_DB(target_db, "so2_max_cumulative", alt_ent_class_target, value)
        if emissionCap["entity_byname"][1] == "NOX" or emissionCap["entity_byname"][1] == "nox":
            target_db = ines_transform.add_item_to_DB(target_db, "nox_max_cumulative", alt_ent_class_target, value)
    for emissionTax in emissionTaxs:
        value = api.from_database(emissionTax["value"], emissionTax["type"])
        alt_ent_class_target = [emissionTax["alternative_name"],(emissionTax["entity_byname"][0],), "group"]
        if emissionTax["entity_byname"][1] == "CO2" or emissionTax["entity_byname"][1] == "co2":
            target_db = ines_transform.add_item_to_DB(target_db, "co2_price", alt_ent_class_target, value)
        if emissionTax["entity_byname"][1] == "SO2" or emissionTax["entity_byname"][1] == "so2":
            target_db = ines_transform.add_item_to_DB(target_db, "so2_price", alt_ent_class_target, value)
        if emissionTax["entity_byname"][1] == "NOX" or emissionTax["entity_byname"][1] == "nox":
            target_db = ines_transform.add_item_to_DB(target_db, "nox_price", alt_ent_class_target, value)
    
    for fomEmission in fomEmissions:
        value = api.from_database(fomEmission["value"], fomEmission["type"])
        alt_ent_class_target = [fomEmission["alternative_name"],(fomEmission["entity_byname"][0],), "unit"]
        if emissionTax["entity_byname"][1] == "CO2" or emissionTax["entity_byname"][1] == "co2":
            target_db = ines_transform.add_item_to_DB(target_db, "fixed_co2_emissions", alt_ent_class_target, value)
    
    for invEmission in invEmissions:
        value = api.from_database(invEmission["value"], invEmission["type"])
        alt_ent_class_target = [invEmission["alternative_name"],(invEmission["entity_byname"][0],), "unit"]
        if emissionTax["entity_byname"][1] == "CO2" or emissionTax["entity_byname"][1] == "co2":
            target_db = ines_transform.add_item_to_DB(target_db, "investment_co2_emissions", alt_ent_class_target, value)
    for param in emissionPriceChanges:
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"], (param["entity_byname"][0],),'set']
        target_db = single_price_change(target_db, t_val__timestamp, value, alt_ent_class_target, 'co2_price', 'co2_price_forecasts')


    for vomEmission in vomEmissions:
        value = api.from_database(vomEmission["value"], vomEmission["type"])
        for gnui in grid__node__unit__io:
            if vomEmission["entity_byname"] == (gnui["entity_byname"][1],gnui["entity_byname"][2]):
                if gnui["entity_byname"][3] == "input":
                    alt_ent_class_target = [vomEmission["alternative_name"],
                                            (vomEmission["entity_byname"][1],vomEmission["entity_byname"][2]), "node__to_unit"]
                    target_db = ines_transform.add_item_to_DB(target_db, "nox_price", alt_ent_class_target, value)
                elif gnui["entity_byname"][3] == "output":
                    alt_ent_class_target = [vomEmission["alternative_name"],
                                            (vomEmission["entity_byname"][2],vomEmission["entity_byname"][1]), "unit__to_node"]
                if vomEmission["entity_byname"][3] == "SO2" or vomEmission["entity_byname"][3] == "so2":
                    target_db = ines_transform.add_item_to_DB(target_db, "so2_emission_rate", alt_ent_class_target, value)
                if vomEmission["entity_byname"][3] == "NOX" or vomEmission["entity_byname"][3] == "nox":
                    target_db = ines_transform.add_item_to_DB(target_db, "nox_emission_rate", alt_ent_class_target, value)    
    
    return target_db

def create_reserves(source_db, target_db, t_val__timestamp):

    restypes = source_db.get_entity_items(entity_class_name='restype')
    group__restypes = source_db.get_entity_items(entity_class_name='group__restype')
    group__restype__uds = source_db.get_entity_items(entity_class_name='group__restype__up_down')
    group__nodes = source_db.get_entity_items(entity_class_name='group__node')
    gnurs = source_db.get_entity_items(entity_class_name='grid__node__unit__restype')
    gnnrs = source_db.get_entity_items(entity_class_name='grid__node__node__restype')
    reserve_activation_durations = source_db.get_parameter_value_items(entity_class_name='group__restype', parameter_definition_name = 'reserve_activation_duration') 
    reserve_lengths = source_db.get_parameter_value_items(entity_class_name='group__restype', parameter_definition_name = 'reserve_length')
    gate_closures = source_db.get_parameter_value_items(entity_class_name='group__restype', parameter_definition_name = 'gate_closure')
    reserve_reactivation_times = source_db.get_parameter_value_items(entity_class_name='group__restype', parameter_definition_name = 'reserve_reactivation_time')
    update_frequencies = source_db.get_parameter_value_items(entity_class_name='group__restype', parameter_definition_name = 'update_frequency')
    reserveDemands = source_db.get_parameter_value_items(entity_class_name='group__restype__up_down', parameter_definition_name = 'reserveDemand')
    LossOfTrans = source_db.get_parameter_value_items(entity_class_name='group__restype__up_down', parameter_definition_name = 'LossOfTrans')
    portion_of_transfer_to_reserve = source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name = 'portion_of_transfer_to_reserve')
    portion_of_infeed_to_reserve = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__restype', parameter_definition_name = 'portion_of_infeed_to_reserve') 
    unit_fail = source_db.get_parameter_value_items(entity_class_name='unit', parameter_definition_name = 'unit_fail')
    reserveReliability = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__restype', parameter_definition_name = 'reserveReliability') 
    reserve_increase_ratio = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__restype', parameter_definition_name = 'reserve_increase_ratio')
    gnur_ups = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__restype', parameter_definition_name = 'up')
    gnur_downs = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__restype', parameter_definition_name = 'down')
    gnnr_ups =  source_db.get_parameter_value_items(entity_class_name='grid__node__node__restype', parameter_definition_name = 'up')
    gnnr_downs = source_db.get_parameter_value_items(entity_class_name='grid__node__node__restype', parameter_definition_name = 'down')
    
    #update_offset
    #offlineReserveCapability
    #portion_of_infeed_to_reserve   the coefficient is not applied, only the contingecy is made

    #create entities first:
    for entity in group__restypes:
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__reserve', entity_byname=entity["entity_byname"]), warn=True)  
    for gnnr in gnnrs:
        link_name  = f'link_{gnnr["entity_byname"][1]}_{gnnr["entity_byname"][2]}'
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='link__node__reserve', entity_byname=(link_name,gnnr["entity_byname"][1],gnnr["entity_byname"][3])), warn=True)     
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='link__node__reserve', entity_byname=(link_name,gnnr["entity_byname"][2],gnnr["entity_byname"][3])), warn=True)     
    for gnur in gnurs:
        entity_byname = (gnur["entity_byname"][2],gnur["entity_byname"][1],gnur["entity_byname"][3])
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='unit__node__reserve', entity_byname= entity_byname), warn=True)

    for restype in restypes:
        #should work on most cases, but might fail if same reserve has different reserve amounts for up and down
        up = False
        down = False
        conting = False
        for entity in group__restype__uds:
            if entity["entity_byname"][1] == restype["entity_byname"][0]:
                if entity["entity_byname"][2] == "up":
                    up = True
                elif entity["entity_byname"][2] == "down":
                    down = True

        for param in portion_of_infeed_to_reserve:
            for unit_param in unit_fail:
                value = api.from_database(param["value"], param["type"])
                if unit_param["entity_byname"][0] == param["entity_byname"][2]:
                    value = api.from_database(param["value"], param["type"])
                    if param["entity_byname"][3] == restype["entity_byname"][0] and value > 0:
                        conting = True
                        alt_ent_class_target = [param["alternative_name"], (param["entity_byname"][2],param["entity_byname"][1],param["entity_byname"][3]), "unit__node__reserve"]
                        target_db = ines_transform.add_item_to_DB(target_db, "contingency_causing", alt_ent_class_target, True)
    
        for pottr in portion_of_transfer_to_reserve:
            gnn1 = False
            gnn2 = False
            for gr in group__restypes:
                if gr["entity_byname"][1] == restype["entity_byname"][0]:
                    for gn in group__nodes:
                        if gr["entity_byname"][0] == gn["entity_byname"][0]:
                            if pottr["entity_byname"][1] == gn["entity_byname"][1]:
                                gnn1 = True
                            if pottr["entity_byname"][2] == gn["entity_byname"][1]:
                                gnn2 = True
            if gnn1 ^ gnn2:
                value = api.from_database(pottr["value"], pottr["type"])
                link_name = 'link_'+param["entity_byname"][1]+"_"+param["entity_byname"][2]
                alt_ent_class_target = [param["alternative_name"], (link_name, param["entity_byname"][1], restype["entity_byname"][0]), "link__node__reserve"]
                target_db = ines_transform.add_item_to_DB(target_db, "reserve_requirement_factor", alt_ent_class_target, value)
                alt_ent_class_target = [param["alternative_name"], (link_name, param["entity_byname"][2], restype["entity_byname"][0]), "link__node__reserve"]
                target_db = ines_transform.add_item_to_DB(target_db, "reserve_requirement_factor", alt_ent_class_target, value) 
                for param in LossOfTrans:
                    if param["entity_byname"][1] == restype["entity_byname"][0]:
                        alt_ent_class_target = [pottr["alternative_name"], (link_name, param["entity_byname"][1], param["entity_byname"][3]), "link__node__reserve"]
                        target_db = ines_transform.add_item_to_DB(target_db, "contingency_causing", alt_ent_class_target, True)
                        alt_ent_class_target = [pottr["alternative_name"], (link_name, param["entity_byname"][2], param["entity_byname"][3]), "link__node__reserve"]
                        target_db = ines_transform.add_item_to_DB(target_db, "contingency_causing", alt_ent_class_target, True)
                        conting = True
        
        if conting:
            reserve_type = "contingency"
        elif up and down:
            reserve_type = "symmetric"
        elif up:
            reserve_type = "upward"
        elif down:
            reserve_type = "downward"
    
        alt_ent_class_target = [settings["alternative"],  restype["entity_byname"], "reserve"]
        if conting or up or down:
            target_db = ines_transform.add_item_to_DB(target_db, "reserve_type", alt_ent_class_target, reserve_type)

    for param in reserve_activation_durations:
        value = api.from_database(param["value"], param["type"])
        duration = api.Duration(str(int(value*60))+"min")
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "set__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "activation_time", alt_ent_class_target, duration)
    for param in gate_closures:
        value = api.from_database(param["value"], param["type"]) * settings["stepLengthInHours"]
        duration = api.Duration(str(int(value*60))+"min")
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "set__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "gate_closure", alt_ent_class_target, duration)
    for param in reserve_reactivation_times:
        value = api.from_database(param["value"], param["type"])
        duration = api.Duration(str(int(value*60))+"min")
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "set__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_reactivation_time", alt_ent_class_target, duration) 
    for param in update_frequencies:
        value = api.from_database(param["value"], param["type"])* settings["stepLengthInHours"]
        duration = api.Duration(str(int(value*60))+"min")
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "set__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "update_frequency", alt_ent_class_target, duration)  
    for param in reserve_lengths:
        value = api.from_database(param["value"], param["type"]) * settings["stepLengthInHours"]
        alt_ent_class_target = [param["alternative_name"], param["entity_byname"], "set__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "duration", alt_ent_class_target, duration) 
    for param in reserveDemands:
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][0], param["entity_byname"][1]), "set__reserve"]
        target_db = pass_timeseries(target_db, 'reserve_requirement', 'reserve_requirement_forecasts', value, alt_ent_class_target, t_val__timestamp)
    
    for param in reserve_increase_ratio:
        value = api.from_database(param["value"], param["type"])
        alt_ent_class_target = [param["alternative_name"],(param["entity_byname"][2], param["entity_byname"][1],param["entity_byname"][3]), "unit__node__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_requirement_factor", alt_ent_class_target, value)  

    for param in gnur_ups:
        value = api.from_database(param["value"], param["type"])
        for rR in reserveReliability:
            if param["entity_byname"] == rR["entity_byname"]:
                value = value * api.from_database(rR["value"], rR["type"])
        alt_ent_class_target = [param["alternative_name"] ,(param["entity_byname"][2], param["entity_byname"][1],param["entity_byname"][3]), "unit__node__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_provision_coefficient", alt_ent_class_target, value)
    
    for param in gnur_downs:
        value = api.from_database(param["value"], param["type"])
        for rR in reserveReliability:
            if param["entity_byname"] == rR["entity_byname"]:
                value = value * api.from_database(rR["value"], rR["type"])
        alt_ent_class_target = [param["alternative_name"] ,(param["entity_byname"][2], param["entity_byname"][1],param["entity_byname"][3]), "unit__node__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_provision_coefficient", alt_ent_class_target, value)

    #Is the reserve actually group or group-restype not just restype
    #are links in reserves only on the border of the group
    for param in gnnr_ups:
        value = api.from_database(param["value"], param["type"])
        link_name  = f'link_{param["entity_byname"][1]}_{param["entity_byname"][2]}'
        alt_ent_class_target = [param["alternative_name"] , (link_name, param["entity_byname"][1], param["entity_byname"][3]), "link__node__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_provision_coefficient", alt_ent_class_target, value)
    
    for param in gnnr_downs:
        value = api.from_database(param["value"], param["type"])
        link_name  = f'link_{param["entity_byname"][1]}_{param["entity_byname"][2]}'
        alt_ent_class_target = [param["alternative_name"] , (link_name, param["entity_byname"][1], param["entity_byname"][3]), "link__node__reserve"]
        target_db = ines_transform.add_item_to_DB(target_db, "reserve_provision_coefficient", alt_ent_class_target, value)
    
    return target_db

def create_unit_node_constraints(source_db, target_db, t_val__timestamp):

    gnuios = source_db.get_entity_items(entity_class_name='grid__node__unit__io')
    constants = source_db.get_parameter_value_items(entity_class_name='unit__constraint', parameter_definition_name = 'constant')
    unit_constraints = source_db.get_entity_items(entity_class_name='unit__constraint')
    coefficients = source_db.get_parameter_value_items(entity_class_name='unit__constraint__node', parameter_definition_name = 'coefficient')
    units = list()

    for unit_constraint in unit_constraints:
        if unit_constraint["entity_byname"][0] not in units:
            units.append(unit_constraint["entity_byname"][0])
            
        uc = unit_constraint["entity_byname"]
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', 
                                                                    entity_byname=(f"{uc[0]}_{uc[1]}",)), warn=True)
        constant_value = 0
        for constant in constants:
            if constant["entity_byname"] == uc:
                constant_value = api.from_database(constant["value"], constant["type"])
        target_db = ines_transform.add_item_to_DB(target_db, "constant", [settings['alternative'],(f"{uc[0]}_{uc[1]}",),'constraint'], constant_value) 

        if uc[1][0:2] == "eq":
            sense = "equal"
        elif uc[1][0:2] == "gt": 
            sense = "greater_than"
        else:
            sense = "less_than"
        target_db = ines_transform.add_item_to_DB(target_db, "sense", [settings['alternative'],(f"{uc[0]}_{uc[1]}",),'constraint'], sense) 

    
    for unit in units:            
        for gnuio in gnuios:
            if gnuio["entity_byname"][2] == unit:
                value_indexes = list()
                value_values = list()
                for coefficient in coefficients:
                    if gnuio["entity_byname"][1] == coefficient["entity_byname"][2] and gnuio["entity_byname"][2] == coefficient["entity_byname"][0]:
                        coefficient_value = api.from_database(coefficient["value"], coefficient["type"])
                        value_indexes.append(f"{coefficient["entity_byname"][0]}_{coefficient["entity_byname"][1]}")
                        value_values.append(coefficient_value)
                        alt = coefficient['alternative_name']
                if len(value_indexes) > 0: 
                    value = api.Map(value_indexes, value_values)
                    if gnuio["entity_byname"][3] == 'input':
                        target_db = ines_transform.add_item_to_DB(target_db, "constraint_flow_coefficient", [alt,(gnuio["entity_byname"][1],gnuio["entity_byname"][2]),'node__to_unit'], value) 
                    elif gnuio["entity_byname"][3] == 'output':
                        target_db = ines_transform.add_item_to_DB(target_db, "constraint_flow_coefficient", [alt,(gnuio["entity_byname"][2],gnuio["entity_byname"][1]),'unit__to_node'], value) 

    return target_db

def create_node_capacities(source_db, target_db, t_val__timestamp):
    
    #backbone does not have node capacities. Instead the node state is limited by upperLimit (MWh, not 0-1).
    #Investments can be made only through grid__node__unit__io with upperLimitCapacityRatio. This increases the upperLimit by summing to it.

    #this should be rewritten
    #some notes:
        #all of unit__node__boundary upperLimit and lowerLimit exist even if they have no parameters


    nodes = source_db.get_entity_items(entity_class_name='node')
    gnbs = source_db.get_entity_items(entity_class_name='grid__node__boundary')
    energyStoredPerUnitOfStates = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'energyStoredPerUnitOfState')
    constants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'constant')
    multipliers = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'multiplier')
    useConstants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useConstant')
    timeseriess = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'timeseries')
    useTimeseriess = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useTimeSeries')
    unit_sizes = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name = 'unitSize')
    unit_capacities = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name = 'capacity')
    unit_counts = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name = 'unitCount')
    uLCRs = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__io', parameter_definition_name = 'upperLimitCapacityRatio')
    
    for node in nodes:
        out = 0
        eSPUOS = 1
        for energyStoredPerUnitOfState in energyStoredPerUnitOfStates:
            if node["entity_byname"][0] == energyStoredPerUnitOfState["entity_byname"][1]:
                eSPUOS = api.from_database(energyStoredPerUnitOfState["value"],energyStoredPerUnitOfState["type"])
        # first search for the capacity from upward limit, the lower limit is also per unit so the capacity is needed first
        capacity = 0
        capacity_found = False
        for gnb in gnbs:
            multi = 1
            for multiplier in multipliers:
                if gnb["entity_byname"] == multiplier["entity_byname"]:
                    multi = api.from_database(multiplier["value"], multiplier["type"])
            if node["entity_byname"][0] == gnb["entity_byname"][1] and gnb["entity_byname"][2] == "upwardLimit":
                if any(gnb["entity_byname"] == useConstant["entity_byname"] for useConstant in useConstants):
                    for constant in constants:
                        if gnb["entity_byname"] == constant["entity_byname"]:
                            capacity = api.from_database(constant["value"], constant["type"]) * multi*eSPUOS
                            alt = constant["alternative_name"]
                elif any(gnb["entity_byname"] == useTimeseries["entity_byname"] for useTimeseries in useTimeseriess):
                    for timeseries in timeseriess:
                        if gnb["entity_byname"] == timeseries["entity_byname"]:
                            out = api.from_database(timeseries["value"], timeseries["type"])
                            alt = timeseries["alternative_name"]
                            #search for the max value in the timeseries, two or one dimensional
                            if isinstance(out.values[0], api.parameter_value.Map):
                                for values_map in out.values:
                                    for val in values_map.values:
                                        val = val * multi *eSPUOS
                                        if val > capacity:
                                            capacity = val
                            else:
                                for val in out.values:
                                    val = val * multi *eSPUOS
                                    if val > capacity:
                                            capacity = val
                capacity_found = True
        for uLCR in uLCRs:
            if node["entity_byname"][0] == uLCR["entity_byname"][1]:
                alt = uLCR["alternative_name"]
                #two ways of representing unit capacities
                for unit_size in unit_sizes:
                    if unit_size["entity_byname"] == uLCR["entity_byname"]:
                        unit_size_value = api.from_database(unit_size["value"],unit_size["type"])
                        ratio = api.from_database(uLCR["value"],uLCR["type"])
                        u_count = 1
                        for unit_count in unit_counts:
                            if unit_count["entity_byname"] == uLCR["entity_byname"]:
                                u_count = api.from_database(unit_count["value"],unit_count["type"])
                        capacity = capacity + ratio * unit_size_value * u_count  * eSPUOS
                        capacity_found = True
                for unit_capacity in unit_capacities:
                    if unit_capacity["entity_byname"] == uLCR["entity_byname"]:
                        unit_capacity_value = api.from_database(unit_capacity["value"],unit_capacity["type"])
                        capacity = capacity + unit_capacity_value * eSPUOS
                        capacity_found = True
        
        if capacity_found:
            node_alt_ent_class_target = [alt, (node["entity_byname"][0],), "node"]
            target_db = ines_transform.add_item_to_DB(target_db, 'storages_existing', node_alt_ent_class_target, 1)
            target_db = ines_transform.add_item_to_DB(target_db, 'storage_capacity', node_alt_ent_class_target, capacity)

        # Add upper_limit timeseries                 
        for gnb in gnbs:
            multi = 1
            for multiplier in multipliers:
                if gnb["entity_byname"] == multiplier["entity_byname"]:
                    multi = api.from_database(multiplier["value"], multiplier["type"])
            if node["entity_byname"][0] == gnb["entity_byname"][1] and gnb["entity_byname"][2] == "upwardLimit":
                if any(gnb["entity_byname"] == useTimeseries["entity_byname"] for useTimeseries in useTimeseriess):
                    for timeseries in timeseriess:
                        if gnb["entity_byname"] == timeseries["entity_byname"]:
                            out = api.from_database(timeseries["value"], timeseries["type"])                    
                            if isinstance(out.values[0], api.parameter_value.Map):
                                for i, values_map in enumerate(out.values):
                                    for j, val in enumerate(values_map.values):
                                        out.values[i].values[j] = val * multi *eSPUOS / capacity
                            else:
                                for i, val in enumerate(out.values):
                                    out.values[i] = val * multi *eSPUOS / capacity
                            if not any(gnb["entity_byname"][1] == uLCR["entity_byname"][1] for uLCR in uLCRs):
                                target_db = pass_timeseries(target_db, 'storage_state_upper_limit', 'storage_state_upper_limit_forecasts', out, node_alt_ent_class_target, t_val__timestamp)
        
        # add the downward limit next
        for gnb in gnbs:
            multi = 1
            for multiplier in multipliers:
                if gnb["entity_byname"] == multiplier["entity_byname"]:
                    multi = api.from_database(multiplier["value"], multiplier["type"])
            if node["entity_byname"][0] == gnb["entity_byname"][1] and (gnb["entity_byname"][2] == "downwardLimit" or gnb["entity_byname"][2] == "reference") :
                if any(gnb["entity_byname"] == useConstant["entity_byname"] for useConstant in useConstants):
                    for constant in constants:
                        if gnb["entity_byname"] == constant["entity_byname"]:
                            value = api.from_database(constant["value"], constant["type"])
                            node_alt_ent_class_target = [constant["alternative_name"], (gnb["entity_byname"][1],), "node"]
                            if gnb["entity_byname"][2] == "downwardLimit":
                                target_db = ines_transform.add_item_to_DB(target_db, 'storage_state_lower_limit',  node_alt_ent_class_target, value * multi * eSPUOS/capacity)
                elif any(gnb["entity_byname"] == useTimeseries["entity_byname"] for useTimeseries in useTimeseriess):
                    for timeseries in timeseriess:
                        if gnb["entity_byname"] == timeseries["entity_byname"]:
                            out = api.from_database(timeseries["value"], timeseries["type"])
                            node_alt_ent_class_target = [timeseries["alternative_name"], (gnb["entity_byname"][1],), "node"]
                            if isinstance(out.values[0], api.parameter_value.Map):
                                for i, values_map in enumerate(out.values):
                                    for j, val in enumerate(values_map.values):
                                        out.values[i].values[j] = val * multi *eSPUOS / capacity
                            else:
                                for i, val in enumerate(out.values):
                                    out.values[i] = val * multi *eSPUOS / capacity
                            if gnb["entity_byname"][2] == "downwardLimit":
                                target_db = pass_timeseries(target_db, 'storage_state_lower_limit', "storage_state_lower_limit_forecasts", out, node_alt_ent_class_target, t_val__timestamp)
                            if gnb["entity_byname"][2] == "reference":
                                target_db = pass_timeseries(target_db, 'storage_state_fix', 'storage_state_fix_forecasts', out, node_alt_ent_class_target, t_val__timestamp)
        for uLCR in uLCRs:
            size_found = False
            capacity_found = False
            ratio = api.from_database(uLCR["value"],uLCR["type"])
            if node["entity_byname"][0] == uLCR["entity_byname"][1]:
                for unit_size in unit_sizes:
                    if unit_size["entity_byname"] == uLCR["entity_byname"]:
                        unit_size_value = api.from_database(unit_size["value"],unit_size["type"])
                        size_found = True
                u_count = 1
                for unit_count in unit_counts:
                    if unit_count["entity_byname"] == uLCR["entity_byname"]:
                        u_count = api.from_database(unit_count["value"],unit_count["type"])
                for unit_capacity in unit_capacities:
                    if unit_capacity["entity_byname"] == uLCR["entity_byname"]:
                        unit_capacity_value = api.from_database(unit_capacity["value"],unit_capacity["type"])
                    capacity_found = True

                constraint_byname = (uLCR["entity_byname"][1]+"_"+uLCR["entity_byname"][2]+"_investment",)
                unit_alt_ent_class_target = [uLCR["alternative_name"], (uLCR["entity_byname"][2],), "unit"]
                node_alt_ent_class_target = [uLCR["alternative_name"], (uLCR["entity_byname"][1],), "node"]
                target_db = ines_transform.add_item_to_DB(target_db, 'storage_investment_method', node_alt_ent_class_target, "no_limits")
                target_db = ines_transform.add_item_to_DB(target_db, "storage_investment_cost", node_alt_ent_class_target, 0)
                
                #constraint
                if size_found:
                    const = 1 - (ratio * unit_size_value * u_count) / capacity
                    unit_coeff = - ratio * unit_size_value / capacity
                elif capacity_found:
                    const = 1 - (ratio * unit_capacity_value) / capacity
                    unit_coeff = - ratio * unit_capacity_value / capacity
                else:
                    print("Warning: grid__node__unit__io with upperLimitCapacityRatio without a capacity defined")
                    continue
                ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', entity_byname=constraint_byname), warn=True)
                target_db = ines_transform.add_item_to_DB(target_db, 'sense', [uLCR["alternative_name"], constraint_byname, "constraint"], 'greater_than')
                target_db = ines_transform.add_item_to_DB(target_db, 'constant', [uLCR["alternative_name"], constraint_byname, "constraint"], const)
                target_db = ines_transform.add_item_to_DB(target_db, 'constraint_storage_count_coefficient', node_alt_ent_class_target, api.Map([constraint_byname[0]],[1]))
                target_db = ines_transform.add_item_to_DB(target_db, 'constraint_unit_count_coefficient', unit_alt_ent_class_target, api.Map([constraint_byname[0]],[unit_coeff]))

    return target_db

def create_bound_state_constraints(source_db, target_db, t_val__timestamp):
    #not tested
    for param in source_db.get_parameter_value_items(entity_class_name='grid__node__node', parameter_definition_name = 'boundStateMaxDiff'):
        value = api.from_database(param["value"], param["type"])
        constraint_map = {
            "less_than": ['state_leq',param["entity_byname"][1], param["entity_byname"][2]],
            "greater_than": ['state_geq', param["entity_byname"][2], param["entity_byname"][1]]
        }
        for name, value in constraint_map:
            target_byname = ("_".join(value[0],value[1],value[2]),)
            ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', 
                                                                    entity_byname=target_byname), warn=True)
            alt_ent_class_constraint = [param["alternative_name"], target_byname, "constraint"]
            target_db = ines_transform.add_item_to_DB(target_db, 'sense', alt_ent_class_constraint, name)
            target_db = ines_transform.add_item_to_DB(target_db, 'constant', alt_ent_class_constraint, value)
            
            node_param = api.Map([target_byname], [1])
            alt_ent_class_node = [param["alternative_name"],(param["entity_byname"][1],), "node"]
            target_db = ines_transform.add_item_to_DB(target_db, 'constraint_storage_state_coefficient', alt_ent_class_node, node_param)
            node_param = api.Map([target_byname], [-1])
            alt_ent_class_node = [param["alternative_name"],(param["entity_byname"][2],), "node"]
            target_db = ines_transform.add_item_to_DB(target_db, 'constraint_storage_state_coefficient', alt_ent_class_node, node_param)
    
    return target_db

def handle_boundaries(source_db, target_db, t_val__timestamp):
  
    gnubs = source_db.get_entity_items(entity_class_name='grid__node__unit__boundary')
    rampCosts = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__boundary', parameter_definition_name = 'rampCost')
    rampLimits = source_db.get_parameter_value_items(entity_class_name='grid__node__unit__boundary', parameter_definition_name = 'rampLimit')
    gnuios = source_db.get_entity_items(entity_class_name='grid__node__unit__io')


    gnbs = source_db.get_entity_items(entity_class_name='grid__node__boundary')
    constants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'constant')
    multipliers = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'multiplier')
    timeseriess = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'timeseries')
    useConstants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useConstant')
    useTimeseriess = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useTimeseries')
    energyStoredPerUnitOfStates = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'energyStoredPerUnitOfState')

    for gnb in gnbs:
        eSPUOS = 1
        for energyStoredPerUnitOfState in energyStoredPerUnitOfStates:
            if energyStoredPerUnitOfState["entity_byname"][1] == gnb["entity_byname"][1]:
                eSPUOS = api.from_database(energyStoredPerUnitOfState["value"],energyStoredPerUnitOfState["type"])
        multi = 1
        for multiplier in multipliers:
            if gnb["entity_byname"] == multiplier["entity_byname"]:
                multi = api.from_database(multiplier["value"], multiplier["type"])
        if any(gnb["entity_byname"] == useConstant["entity_byname"] for useConstant in useConstants):
            for constant in constants:
                if gnb["entity_byname"] == constant["entity_byname"]:
                    out = api.from_database(constant["value"], constant["type"]) * multi * eSPUOS
                    alt = constant["alternative_name"]
        elif any(gnb["entity_byname"] == useTimeseries["entity_byname"] for useTimeseries in useTimeseriess):
            for timeseries in timeseriess:
                #check if works with stochastic timeseries if they exist here
                if gnb["entity_byname"] == timeseries["entity_byname"]:
                    out = api.from_database(timeseries["value"], timeseries["type"])
                    if isinstance(out.values[0], api.parameter_value.Map):
                        for i, values_map in enumerate(out.values):
                            for j, val in enumerate(values_map.values):
                                out.values[i].values[j] = val * multi * eSPUOS
                    else:
                        for i, val in enumerate(out.values):
                            out.values[i] = val * multi * eSPUOS
                    alt = timeseries["alternative_name"]

        if gnb["entity_byname"][2] == "maxSpill":
            target_db = pass_timeseries(target_db, 'spill_upper_limit', None, out, [alt, (gnb["entity_byname"][1],) ,"node"], t_val__timestamp)
        elif gnb["entity_byname"][2] == "minSpill":
            target_db = pass_timeseries(target_db, 'spill_lower_limit', None, out, [alt, (gnb["entity_byname"][1],) ,"node"], t_val__timestamp)
        elif gnb["entity_byname"][2] == "balancePenalty":
            target_db = pass_timeseries(target_db, 'penalty_downward', None, out, [alt, (gnb["entity_byname"][1],) ,"node"], t_val__timestamp)
            target_db = pass_timeseries(target_db, 'penalty_upward', None, out, [alt, (gnb["entity_byname"][1],) ,"node"], t_val__timestamp)

    for gnub in gnubs:
        for gnuio in gnuios:
            if gnuio["entity_byname"][1] == gnub["entity_byname"][1] and gnuio["entity_byname"][2] == gnub["entity_byname"][2]:
                if gnuio["entity_byname"][3] == 'input':
                    ent = (gnuio["entity_byname"][1],gnuio["entity_byname"][2])
                    cla = "node__to_unit"
                elif gnuio["entity_byname"][3] == 'output':
                    ent = (gnuio["entity_byname"][2],gnuio["entity_byname"][1])
                    cla = "unit__to_node"

                for rampCost in rampCosts:
                    if gnub["entity_byname"] == rampCost["entity_byname"]:
                        value = api.from_database(rampCost["value"], rampCost["type"])
                        target_db = ines_transform.add_item_to_DB(target_db, 'ramp_cost', [rampCost["alternative_name"], ent, cla], value)

                for rampLimit in rampLimits:
                    if gnub["entity_byname"] == rampLimit["entity_byname"]:
                        value = api.from_database(rampLimit["value"], rampLimit["type"])
                        target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_down', [rampLimit["alternative_name"], ent, cla], value)
                        target_db = ines_transform.add_item_to_DB(target_db, 'ramp_limit_up', [rampLimit["alternative_name"], ent, cla], value)
        #downwardSlack01 - downwardSlack20
        #upwardSlack01 - upwardSlack20
    return target_db

def create_group_constraints(source_db, target_db):

    unit__groups =  source_db.get_entity_items(entity_class_name='unit__group')
    constrainedOnlineMultipliers = source_db.get_parameter_value_items(entity_class_name='unit__group', parameter_definition_name = 'constrainedOnlineMultiplier')
    constrainedCapMultipliers = source_db.get_parameter_value_items(entity_class_name='unit__group', parameter_definition_name = 'constrainedCapMultiplier')
    for param in source_db.get_parameter_value_items(entity_class_name='group', parameter_definition_name = 'constrainedOnlineTotalMax'):
        #create constraint
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', 
                                                                    entity_byname=('constraint_online_'+ param["entity_byname"][0],)), warn=True)
        value =  api.from_database(param["value"], param["type"])
        target_db = ines_transform.add_item_to_DB(target_db, "constant", [param['alternative_name'],('constraint_online_'+ param["entity_byname"][0],),'constraint'], value) 
        target_db = ines_transform.add_item_to_DB(target_db, "sense", [param['alternative_name'],('constraint_online_'+ param["entity_byname"][0],),'constraint'], "less_than") 

        #add units to the constraint
        for unit__group in unit__groups:
            if param["entity_byname"][0] == unit__group["entity_byname"][1]:
                constant_value = None
                for constant in constrainedOnlineMultipliers:
                    if constant["entity_byname"] == unit__group["entity_byname"]:
                        constant_value = api.from_database(constant["value"], constant["type"])
                if constant_value:
                    value = api.Map(['constraint_online_'+ param["entity_byname"][0]],[constant_value])
                    target_db = ines_transform.add_item_to_DB(target_db, "constraint_online_coefficient", [param['alternative_name'],(unit__group["entity_byname"][0],),'unit'], value) 

    for param in source_db.get_parameter_value_items(entity_class_name='group', parameter_definition_name = 'constrainedCapTotalMax'):
        #check if any multipliers
        value =  api.from_database(param["value"], param["type"])
        if not any(param["entity_byname"][0] == multi["entity_byname"][1] for multi in constrainedOnlineMultipliers):
            target_db = ines_transform.add_item_to_DB(target_db, "invest_max_total", [param['alternative_name'],param["entity_byname"],'set'], value) 
        else:
            #use constraint entities
            ines_transform.assert_success(target_db.add_entity_item(entity_class_name='constraint', 
                                                                        entity_byname=('constraint_cap_'+ param["entity_byname"][0],)), warn=True)
            target_db = ines_transform.add_item_to_DB(target_db, "constant", [param['alternative_name'],('constraint_cap_'+ param["entity_byname"][0],),'constraint'], value) 
            target_db = ines_transform.add_item_to_DB(target_db, "sense", [param['alternative_name'],('constraint_cap_'+ param["entity_byname"][0],),'constraint'], "less_than") 

            #add units to the constraint
            for unit__group in unit__groups:
                if param["entity_byname"][0] == unit__group["entity_byname"][1]:
                    constant_value = None
                    for constant in constrainedCapMultipliers:
                        if constant["entity_byname"] == unit__group["entity_byname"]:
                            constant_value = api.from_database(constant["value"], constant["type"])
                    if constant_value:
                        value = api.Map(['constraint_cap_'+ param["entity_byname"][0]],[constant_value])
                        target_db = ines_transform.add_item_to_DB(target_db, "constraint_cap_coefficient", [param['alternative_name'],(unit__group["entity_byname"][0],),'unit'], value) 

    return target_db

def add_node_types(source_db, target_db):

    grid_node_boudaries =  source_db.get_entity_items(entity_class_name='grid__node__boundary')
    use_prices = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'usePrice')
    prices = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'price')
    price_changes = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'priceChange')
    node_balances = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'nodeBalance')
    use_constants = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useConstant')
    use_timeseries = source_db.get_parameter_value_items(entity_class_name='grid__node__boundary', parameter_definition_name = 'useTimeSeries')
    

    for node in source_db.get_entity_items(entity_class_name='node'):
        node_type = None
        found = False
        alt = settings["alternative"]
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
            for use_price in use_prices:
                if use_price["entity_byname"][1] == node["name"]:
                    if api.from_database(use_price["value"],use_price["type"]) != 0.0:
                        node_type = "commodity"
                        alt = use_price["alternative_name"]
                        found = True
                        break
        if not found:
            for grid_node_boudary in grid_node_boudaries:
                if ((grid_node_boudary["entity_byname"][2] == 'upwardLimit' or 
                        grid_node_boudary["entity_byname"][2] == 'downwardLimit') and
                        grid_node_boudary["entity_byname"][1] == node["name"]):
                    for use_constant in use_constants:
                        if use_constant["entity_name"] == grid_node_boudary["name"] and api.from_database(use_constant["value"],use_constant["type"]) != 0.0:
                            node_type = "storage" 
                            alt =  use_constant["alternative_name"]
                            found = True
                            break
                    if not found:
                        for use_timeserie in use_timeseries:
                            if use_timeserie["entity_name"] == grid_node_boudary["name"] and api.from_database(use_timeserie["value"],use_timeserie["type"]) != 0.0:
                                node_type = "storage" 
                                alt = use_timeserie["alternative_name"]
                                found = True
                                break  
        if not found:
            for node_balance in node_balances:
                if node_balance["entity_byname"][1]== node["name"] and api.from_database(node_balance["value"],node_balance["type"]) != 0.0:
                    node_type = "balance"
                    alt = node_balance["alternative_name"]
                    break

        target_db = ines_transform.add_item_to_DB(target_db, 'node_type', [alt, node["entity_byname"], "node"], node_type)
    
    return target_db

def create_simple_timeseries(source_db, target_db, t_val__timestamp):

    influxes = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'influx')
    c_price = source_db.get_parameter_value_items(entity_class_name='node', parameter_definition_name = 'price')
    storage_value = source_db.get_parameter_value_items(entity_class_name='grid__node', parameter_definition_name = 'storageValue')

    for influx in influxes:
        value = api.from_database(influx["value"], influx["type"])
        target_db = pass_timeseries(target_db, 'flow_profile', 'flow_profile_forecasts', value, [influx["alternative_name"], (influx["entity_byname"][1],), "node"], t_val__timestamp)
        target_db = ines_transform.add_item_to_DB(target_db, "flow_scaling_method", [influx["alternative_name"], (influx["entity_byname"][1],), "node"], "use_profile_directly")
    
    for param in c_price:
        value =  api.from_database(param["value"], param["type"])
        target_db = pass_timeseries(target_db, 'commodity_price', 'commodity_price_forecasts', value, [param["alternative_name"], (param["entity_byname"][0],), "node"], t_val__timestamp)
    
    for param in storage_value:
        value =  api.from_database(param["value"], param["type"])
        target_db = pass_timeseries(target_db, 'storage_state_value', 'storage_state_value_forecasts', value, [param["alternative_name"], (param["entity_byname"][0],), "node"], t_val__timestamp)

    return target_db

def create_sets_from_grids(source_db, target_db):
    grids = source_db.get_entity_items(entity_class_name='grid')
    grid__nodes = source_db.get_entity_items(entity_class_name='grid__node')
    grid__node__nodes = source_db.get_entity_items(entity_class_name='grid__node__node')
    grid__node__unit__ios = source_db.get_entity_items(entity_class_name='grid__node__unit__io')

    for grid in grids:
        ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set',entity_byname=("grid_"+grid["name"],)), warn=True)
        for grid__node in grid__nodes:
            if grid__node["entity_byname"][0] == grid["entity_byname"][0]:
                ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__node',
                                                                        entity_byname=("grid_"+grid__node["entity_byname"][0], grid__node["entity_byname"][1])), warn=True)
        for grid__node__node in grid__node__nodes:
            if grid__node__node["entity_byname"][0] == grid["entity_byname"][0]:
                link_name = f'link_{grid__node__node["entity_byname"][1]}_{grid__node__node["entity_byname"][2]}'
                ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__link',
                                                                        entity_byname=("grid_"+ grid__node__node["entity_byname"][0], link_name)), warn=True)

        for gnuio in grid__node__unit__ios:
            if gnuio["entity_byname"][0] == grid["entity_byname"][0]:
                if gnuio["entity_byname"][3] == 'input':
                    ent = ("grid_"+ gnuio["entity_byname"][0], gnuio["entity_byname"][1],gnuio["entity_byname"][2])
                elif gnuio["entity_byname"][3] == 'output':
                    ent = ("grid_"+ gnuio["entity_byname"][0],gnuio["entity_byname"][2],gnuio["entity_byname"][1])
                ines_transform.assert_success(target_db.add_entity_item(entity_class_name='set__unit_flow',
                                                                        entity_byname=ent), warn=True)

    return target_db

def add_entity_alternative_items(target_db):
    for entity_class in target_db.get_entity_class_items():
        for entity in target_db.get_entity_items(entity_class_name=entity_class["name"]):
            item = target_db.get_entity_alternative_item(
                entity_class_name=entity_class["name"],
                entity_byname=entity["entity_byname"],
                alternative_name=settings["alternative"],
            )
            if not item:
                ines_transform.assert_success(target_db.add_update_entity_alternative_item(
                    entity_class_name=entity_class["name"],
                    entity_byname=entity["entity_byname"],
                    alternative_name=settings["alternative"],
                    active=True,
                ))
    return target_db


#not jet functional
def inflow_timeseries_from_csv(target_db, input_folder, t_val__timestamp):

    input_folder = Path(f"./{input_folder}")
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            df = pd.read_csv(os.path.join(input_folder,filename),header = 0)
            columns = df.columns[3:]
            dff =  df['f']
            forecasts = dff.drop_duplicates().values.tolist()
            for node in columns.values:
                inner_series = []
                indexes = []
                for forecast in forecasts:
                    entity_byname = (node,)
                    df2 = pd.DataFrame(data=df)
                    df2.filter(items=[forecast], axis=0)
                    t_vals = df2['t'].values.tolist()
                    values = df2[node].values.tolist()
                    timestamps = [stamp for stamp in t_val__timestamp.values()]
                    out_vals = [values[t_vals.index(t_val)] for t_val in t_val__timestamp.keys() if t_val in t_vals]
                    indexes.append(forecast)
                    inner_series.append(api.TimeSeriesVariableResolution(timestamps, out_vals, ignore_year = False, repeat=False, index_name="time step"))
                value = api.Map(indexes,inner_series)
                target_db = ines_transform.add_item_to_DB(target_db, 'flow_profile_forecasts', [settings["alternative"], entity_byname, "node"], value)
    return target_db

###############
#HELPPER FUNCTIONS
###############

def make_duration(param):
    value = api.from_database(param["value"], param["type"])
    duration = api.Duration(str(value)+"h")
    return duration

def calculate_investment_cost(source_db, target_db, alt_ent_class_source, alt_ent_class_target, storage = False):
    annuity = ines_transform.get_parameter_from_DB(source_db, "annuity", alt_ent_class_source)
    lifetime = settings["default_lifetime"]
    r = settings["default_interest_rate"]
    investment_cost = annuity /(r /(1 - (1 / (1+r))**lifetime))

    if storage:
        target_db = ines_transform.add_item_to_DB(target_db, "storage_investment_cost", alt_ent_class_target, investment_cost)
    else: 
        target_db = ines_transform.add_item_to_DB(target_db, "investment_cost", alt_ent_class_target, investment_cost)
    return target_db

def pass_timeseries(target_db, target_name, target_name_stoch, value, alt_ent_class_target, t_val__timestamp, stochastic_group = None):
    if isinstance(value,float):
        target_db = ines_transform.add_item_to_DB(target_db, target_name, alt_ent_class_target, value)
    else:
        if len(value.values) > 1:
            if target_name_stoch:
                print(f'processing stochastic timeseries:{target_name_stoch}')
                inner_series = []
                forecasts = []
                timestamps = [stamp for stamp in t_val__timestamp.values()]
                for i, val_map in enumerate(value.values):
                    if value.indexes[i] not in settings["p_mfProbability"]:
                        continue
                    t_vals = list(val_map.indexes)
                    out_vals = [val_map.values[t_vals.index(t_val)] if t_val in t_vals else 0.0 for t_val in t_val__timestamp.keys()]
                    timeseries = api.TimeSeriesVariableResolution(timestamps, out_vals, ignore_year = False, repeat=False, index_name="time step")
                    if value.indexes[i] == settings["mf_realization"]:
                        target_db = ines_transform.add_item_to_DB(target_db, target_name, alt_ent_class_target, timeseries, value_type="map")
                    else:
                        inner_series.append(timeseries)
                        forecasts.append(value.indexes[i])
                value.values = inner_series
                value.indexes = forecasts
                target_db = ines_transform.add_item_to_DB(target_db, target_name_stoch, alt_ent_class_target, value)

                #add entity to the stochastic set
                if not stochastic_group:
                    out_byname = ("stochastics",) + alt_ent_class_target[1]
                else:
                    out_byname = (stochastic_group,)  + alt_ent_class_target[1]
                if alt_ent_class_target[2] == "unit__to_node" or alt_ent_class_target[2] == "node__to_unit":
                    out_class = "set__unit_flow"
                else:
                    out_class = f'set__{alt_ent_class_target[2]}'
                target_db.add_entity_item(entity_class_name=out_class, entity_byname=out_byname)
        else: #strip the f00 if it is the only one
            timestamps = [str(stamp) for stamp in t_val__timestamp.values()]
            for i, val_map in enumerate(value.values):
                t_vals = list(val_map.indexes)
                out_vals = [val_map.values[t_vals.index(t_val)] if t_val in t_vals else 0.0 for t_val in t_val__timestamp.keys()]
            time_series = api.TimeSeriesVariableResolution(timestamps, out_vals, ignore_year = False, repeat=False, index_name="time step")
            target_db = ines_transform.add_item_to_DB(target_db, target_name, alt_ent_class_target, time_series)

    return target_db

def single_price_change(target_db, t_val__timestamp, source_value, alt_ent_class_target, target_param, target_name_stoch, stochastic_group = None):
    #stoch
    if isinstance(source_value.values[0], api.parameter_value.Map):
        if len(source_value.indexes) > 1:
            if target_name_stoch:
                inner_price_timeseries = []
                forecasts = []
                for forecast, values_map in source_value.items():
                    values = list()
                    price = 0
                    priceChange_dict = values_map.to_dict()
                    for i in t_val__timestamp.keys():
                        for step__price in priceChange_dict["data"]:
                            if i == step__price[0]:
                                price = step__price[1]
                                break
                        values.append(price)
                    if forecast == 'f00':
                        for i in values.indexes:
                            timestamps.append(str(t_val__timestamp[i]))
                        target_db = ines_transform.add_item_to_DB(target_db, target_param, alt_ent_class_target, price_timeseries, value_type="map")
                    else:
                        inner_price_timeseries.append(api.Map(t_val__timestamp.values(), values))
                        forecasts.append(forecast)
                price_timeseries = api.Map(forecasts, inner_price_timeseries)
                target_db = ines_transform.add_item_to_DB(target_db, target_name_stoch, alt_ent_class_target, price_timeseries, value_type="map")
                if not stochastic_group:
                    out_byname = ("stochastics",) + alt_ent_class_target[1]
                else:
                    out_byname = (stochastic_group,) + alt_ent_class_target[1]
                if alt_ent_class_target[2] == "unit__to_node" or alt_ent_class_target[2] == "node__to_unit":
                    out_class = "set__unit_flow"
                else:
                    out_class = f'set__{alt_ent_class_target[2]}'
                target_db.add_entity_item(entity_class_name=out_class, entity_byname=out_byname)
        else: #remove f00 
            for forecast, values_map in source_value.items():
                values = list()
                price = 0
                priceChange_dict = values_map.to_dict()
                for i in t_val__timestamp.keys():
                    for step__price in priceChange_dict["data"]:
                        if i == step__price[0]:
                            price = step__price[1]
                            break
                    values.append(price)
            timestamps = []
            for i in values.indexes:
                timestamps.append(str(t_val__timestamp[i]))
            time_series = api.TimeSeriesVariableResolution(timestamps, values.values[0], ignore_year = False, repeat=False, index_name="time step")
            target_db = ines_transform.add_item_to_DB(target_db, target_param, alt_ent_class_target, time_series)
    #single value
    elif len(source_value.values) == 1:
        price_timeseries = float(source_value.values[0])
        target_db = ines_transform.add_item_to_DB(target_db, target_param, alt_ent_class_target, price_timeseries)
    #timeseries
    else:   #create a timeseries
        values = list()
        price = 0
        priceChange_dict = source_value.to_dict()
        for i in t_val__timestamp.keys():
            for step__price in priceChange_dict["data"]:
                if i == step__price[0]:
                    price = step__price[1]
            values.append(price)
        time_series = api.TimeSeriesVariableResolution(t_val__timestamp.values(), values, ignore_year = False, repeat=False, index_name="time step")
        target_db = ines_transform.add_item_to_DB(target_db, target_param, alt_ent_class_target, time_series)
    
    return target_db

if __name__ == "__main__":
    developer_mode = False
    if developer_mode:
        # save entities to yaml file
        save_folder = os.path.dirname(__file__)
        conversion_configuration(conversions = [save_folder+'/backbone_to_ines_entities.yaml', save_folder+'/backbone_to_ines_parameters.yaml', save_folder+'/backbone_to_ines_parameter_methods.yaml',
                                             save_folder+'/backbone_to_ines_parameters_to_relationships.yaml'], overwrite=True)
    else:   
        url_db_in = sys.argv[1]
        url_db_out = sys.argv[2]
        settings_path = './backbone_to_ines_settings.yaml'

        # open yaml files
        entities_to_copy,parameter_transforms,parameter_methods, parameters_to_relationships, parameters_to_parameters = conversion_configuration()
        with open(settings_path,'r') as file:
            settings = yaml.safe_load(file) 
        cannot_convert = []

        main()