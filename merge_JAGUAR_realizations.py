#! /usr/bin/env python

from astropy.io import fits
import os
import sys
import argparse
from glob import glob
from collections import OrderedDict
import numpy as np
import json

_SEP = '_'
_FILE_PREFIX = 'JADES'
_SF_PREFIX = 'SF'
_Q_PREFIX = 'Q'
_VERSION_PREFIX = 'v'
_REALIZATION_PREFIX = 'r'
_MOCK = 'mock'
_SUFFIX = '.fits.gz'

if __name__ == '__main__':

    if sys.version_info[0] < 3:
        raise Exception("This script requires using Python version >= 3")

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--folder',
        help="Folder containing the different JAGUAR realizations",
        type=str,
        dest="folder",
        default=os.getcwd()
        )

    parser.add_argument(
        '--version',
        help="Version of the JAGUAR mock",
        type=str,
        dest="version",
        default="1.2"
        )

    parser.add_argument(
        '--JSON-selection',
        help="JSON file containing galaxy selection criteria",
        type=str,
        dest="json_selection"
        )

    args = parser.parse_args()    

    selection_criteria = None
    if args.json_selection is not None:
        with open(args.json_selection) as f:
            selection_criteria = json.load(f, object_pairs_hook=OrderedDict)

    realizations = OrderedDict()
    realizations["folder"] = list()
    realizations["number"] = list()

    folders = glob(os.path.join(args.folder, _REALIZATION_PREFIX + '*/'))
    for _folder in folders:
        realization_folder = os.path.basename(os.path.normpath(_folder))
        try:
            realization_number = int(realization_folder[1:])
        except:
            continue

        realizations["folder"].append(_folder)
        realizations["number"].append(realization_number)
    
    for key, value in realizations.items():
        realizations[key] = np.array(value)

    sorted_indices = np.argsort(realizations["number"])

    for key, value in realizations.items():
        realizations[key] = np.array(value[sorted_indices])

    # First, just read all files to compute the total number of columns
    tot_rows = 0
    for _folder, _number in zip(realizations["folder"], realizations["number"]):
        print("Realization: ", _number)
        for _gal_type in [_SF_PREFIX, _Q_PREFIX]:
            _file_name = _FILE_PREFIX + _SEP + _gal_type + _SEP + _MOCK + _SEP \
                    + _REALIZATION_PREFIX + str(_number) + _SEP + _VERSION_PREFIX + \
                    args.version + _SUFFIX
            _file_name = os.path.join(_folder, _file_name)

            with fits.open(_file_name) as f:
                _table = f[1]
                _n_rows = len(_table.data.field(0))

                mask = np.ones(_n_rows, dtype=bool)
                if selection_criteria is not None:
                    for criterion in selection_criteria:
                        _col = _table.data[criterion["column"]]
                        _oper = criterion["operation"]
                        mask[~eval('_col'+_oper)] = False

                tot_rows += np.sum(mask)

    # New you can actually merge all the tables !
    new_hdu = None
    first_row = 0
    for _folder, _number in zip(realizations["folder"], realizations["number"]):
        print("Realization: ", _number)
        for _gal_type in [_SF_PREFIX, _Q_PREFIX]:
            _file_name = _FILE_PREFIX + _SEP + _gal_type + _SEP + _MOCK + _SEP \
                    + _REALIZATION_PREFIX + str(_number) + _SEP + _VERSION_PREFIX + \
                    args.version + _SUFFIX
            _file_name = os.path.join(_folder, _file_name)

            with fits.open(_file_name) as f:
                _table = f[1]
                _n_rows = len(_table.data.field(0))

                mask = np.ones(_n_rows, dtype=bool)
                if selection_criteria is not None:
                    for criterion in selection_criteria:
                        _col = _table.data[criterion["column"]]
                        _oper = criterion["operation"]
                        mask[~eval('_col'+_oper)] = False

                _n_mask_rows = np.sum(mask)
                orig_cols = _table.columns

                _IDs = np.array(_table.data["ID"], dtype=str)
                _IDs = np.char.add((_REALIZATION_PREFIX + str(_number) + _SEP + _gal_type
                        + _SEP), _IDs)

                orig_cols.del_col("ID")

                str_len = len(max(_IDs, key=len)) + 2
                col_1 = fits.Column(name='ID', format=str(str_len)+'A', array=_IDs)

                if _gal_type == _SF_PREFIX:
                    col_2 = fits.Column(name='SF', format='L', array=np.ones(_n_rows, dtype=bool))
                    col_3 = fits.Column(name='Q', format='L', array=np.zeros(_n_rows, dtype=bool))
                elif _gal_type == _Q_PREFIX: 
                    col_2 = fits.Column(name='SF', format='L', array=np.zeros(_n_rows, dtype=bool))
                    col_3 = fits.Column(name='Q', format='L', array=np.ones(_n_rows, dtype=bool))

                new_cols = fits.ColDefs([col_1, col_2, col_3])
                all_cols = (new_cols + orig_cols)

                new_table = fits.BinTableHDU.from_columns(all_cols)

                if new_hdu is None:
                    new_hdu = fits.BinTableHDU.from_columns(all_cols, nrows=tot_rows)

                last_row = first_row + _n_mask_rows
                for colname in new_hdu.columns.names:
                    if colname in new_table.columns.names:
                        new_hdu.data[colname][first_row:last_row] = new_table.data[colname][mask]

                first_row += _n_mask_rows

    _file_name = _FILE_PREFIX + _SEP + _MOCK + _SEP \
            + _REALIZATION_PREFIX + str(realizations["number"][0]) + _SEP + "to" \
            + _SEP + _REALIZATION_PREFIX + str(realizations["number"][-1]) \
            + _SEP + _VERSION_PREFIX + args.version + _SUFFIX

    new_hdu.writeto(_file_name, overwrite=True)

#with fits.open(fits_table_filename) as hdul1:
#    with fits.open(fits_table_filename) as hdul2:
#        nrows1 = hdul1[1].data.shape[0]
#        nrows2 = hdul2[1].data.shape[0]
#        nrows = nrows1 + nrows2
#        hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
#        for colname in hdul1[1].columns.names:
#            hdu.data[colname][nrows1:] = hdul2[1].data[colname]


