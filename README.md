# ines-backbone
Translation between the ines specification and the Backbone structure

## Status

Functional, but not completely tested.
Uses ines-spec 1.1.1

## use
Get spine-toolbox "https://github.com/spine-tools/Spine-Toolbox", follow its install instructions. The transformation uses spine_db_api.

Get ines-tools "https://github.com/ines-tools/ines-tools", and install it by "pip install ." in its folder.

Get empty ines-spec database from "https://github.com/ines-tools/ines-spec"

### Prerequisite

Transform the backbone gdx file to spine database using the workflow that comes with the backbone.
Add time settings to backbone_to_ines_settings.yaml this includes the information that is in scheduleInits.gms.

### Option 1. Use directly from command line

Make sure that the environment is on where the spine-toolbox installed.

The command to transform is:

```
python backbone_to_ines.py path_to_input_db path_to_output_db
```

where, input_db is the backbone database and output_db is database of the ines-spec.

### Option 2. Use the toolbox project

Open folder backbone_to_ines_workflow as a spine project. 
It contains 3 elements:
+ 1. database for backbone data in spine format
    + main file: Set the path to the backbone_db
+ 2. backbone_to_ines_transform: tool for conversion from backbone to ines 
    + main file = backbone_to_ines.py
    + additional file = backbone_to_ines_settings.yaml
+ 3. database for pypsa data in ines format
    + main file: set the path to the ines_db

## On timeseries

All parameter data needs to be in the gdx file. If timeseries are not, as with northern european model, they can be added to it by running backbone with --debug. Adding timeseries directly to the ines transformation is in the future development plan.

## Notes about the transformation

- The timesteps are in timestamp form in INES. Backbone timesteps are transformed by setting the timestamp of the first timestep in backbone_to_ines_settings.yaml. 

- If timeseries only has one forecast f00, the forecast dimension is removed.

- Nodes have capacities in INES. The backbone format is changed to this using the maximum of UpwardLimit
If investments are allowed with upperLimitCapacityRatio, the node investments are still bound to the unit in INES.

## Information lost in the transformation

### General model parameters

- t_invest: INES assumes investments only at the starts of periods
- discountFactor: INES uses discount rate and time to calculate discount factors. "Freely" given factors cannot be represented explicitly
- realized forecast: waiting for implementation
- Forecast weigths: waiting for implementation   

### Boundary

- upwardSlack and downwardSlack with slackCost

### Reserves

- reservePrice
- restypeReleasedForRealization
- offlineReserveCapability
- updateOffset

### Unit

- BecomeAvailable, BecomeUnavailable
- BoundSamples
- rampSpeedFromMinLoad
- rampSpeedToMinLoad

### Inertia

- unitSizeMVA
- ROCOF
- defaultFrequency
- restype_inertia
- Dynamic inertia as a whole

### Initial system rate

- initialGeneration
- initialInertia
- initialOnlineStatus

### Emissions

- invEmissionFactor
- Only c02 emissions are added to a commodity node, nox and so2 can exist as operational emission

### Misc

- availabilityCapacityMargin
- energyShareMax
- energyShareMin


### z-dimension

### Timeseries

The following parameters are only transformed if they are constants, not timeseries:

- startupCost
- constraint constant
- constraint coefficient

The following parameters will not transform stochastic timeseries:

- MaxSpill
- MinSpill
- penaltyUpward
- penaltyDownward