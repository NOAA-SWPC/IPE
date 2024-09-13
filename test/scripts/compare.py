#!/usr/bin/env python
from netCDF4 import Dataset
from sys import argv
from os.path import exists

def nc_open(f):
    if not exists(f):
        print(f,'not found')
        raise
    return Dataset(f)

def main():
    if len(argv) != 3:
        print('syntax: python compare.py file1 file2')
        return
    a, b = (nc_open(f) for f in argv[1:])
    print('{:<28} {}'.format('variable_name', 'mean_field_difference'))
    for key in a.variables.keys():
        var = a.variables[key][:]
        diff = abs((b.variables[key][:] - var)/var)
        print('{:<28} {:+0.2e}'.format(key, diff.mean()))


if __name__ == '__main__':
    try:
        main()
    except:
        pass
