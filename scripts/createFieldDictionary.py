#!/bin/python

import sys
import yaml
import datetime

if __name__ == '__main__':

    filename = sys.argv[1] if len(sys.argv) > 1  else  '/Volumes/Kea/devel/MOSSCO/code/external/fabm/code/util/standard_variables/variables.yaml'

    if not filename: sys.exit(1)

    with open(filename, 'r') as fid:
        yml = yaml.safe_load(fid)

    entries = []
    for key, value in yml.items():
        d=[{'standard_name': item['name'], 'canonical_units': item['units']}  for i, item in enumerate(value)]
        entries.extend(d)

    entries.append({'standard_name': 'x_velocity_at_10m_above_sea_surface', 'canonical_units': 'm s-1'})
    entries.append({'standard_name': 'y_velocity_at_10m_above_sea_surface', 'canonical_units': 'm s-1'})

    fieldDict={
    'version_number': 0.1,
    'institution': 'Helmholtz-Zentrum Geesthacht Zentrum für Material- und Küstenforschung',
    'source': 'automatically generated from FABM standard variables, with enhancements from MOSSCO',
    'contact': 'Carsten Lemmen <carsten.lemmen@hereon.de>',
    'last_modified': datetime.datetime.now().isoformat(),
    'description': 'Community-based dictionary for shared coupling fields',
    'entries': entries}

# We could also use aliases:
#  - alias: p
#    standard_name: air_pressure
#


    with open('field_dictionary.yaml', 'w') as fid:
        yaml.dump({'field_dictionary': fieldDict}, stream=fid)
