#!/bin/env/python

# Converts a param.nml to global.nml which contains only few variables that
# are needed by the coupled system to be aligned with param.nml

import sys

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

try:
    import f90nml
except:
    print('You need to install the f90nml package (e.g. pip install f90nml)')
    sys.exit(1)

try:
    import yaml
    represent_dict_order = (lambda self, data:
                            self.represent_mapping('tag:yaml.org,2002:map',
                                                   data.items()))
    yaml.add_representer(OrderedDict, represent_dict_order)

except:
    print('You need to install the yaml package (e.g. pip install yaml)')
    sys.exit(1)

if __name__ == '__main__':

    filename = sys.argv[1] if len(sys.argv) > 1 else './param.nml'

    with open(filename) as fid:
        nml = f90nml.read(fid)

    allDict = nml.todict(complex_tuple=True)

    with open(filename.split('.nml')[0] + '.yaml','w') as fid:
        yaml.dump(allDict, fid, default_flow_style=False)

    selectList=[]

    globalDict = {k: v for k, v in allDict['opt'].items() if k.startswith('start_')}
    globalDict['rnday'] = allDict['core']['rnday']

    with open(filename.split('param.nml')[0] + 'global.yaml','w') as fid:
        yaml.dump(globalDict, fid, default_flow_style=False)

    with open(filename.split('param.nml')[0] + 'global.rc','w') as fid:
        yaml.dump(globalDict, fid, default_flow_style=False)

    nml = f90nml.namelist.Namelist({'global': globalDict})

    with open(filename.split('param.nml')[0] + 'global.nml','w') as fid:
        nml.write(fid, force=True, sort=True)
