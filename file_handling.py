"""
Contains input and output functions specific to this project.
"""

import configparser as cp

import pandas as pd

import clt
import gs_opt


# For config files (given as configparser.ConfigParser instances):

def constr_from_config(config: cp.ConfigParser) -> clt.LaminateConstraints:
    return clt.LaminateConstraints(ply_thickness=float(config['PREPREG']['t']),
                                   balanced=config['CONSTR']['balanced'] == 'True',
                                   symmetric=config['CONSTR']['symmetric'] == 'True',
                                   min_ply_share=float(config['CONSTR']['minimum ply share']),
                                   ply_group=[int(ply_angle) for ply_angle in config['CONSTR']['ply group'].split()],
                                   allow_zero_thickness_plies=config['CONSTR']['allow zero-thickness plies'] == 'True',
                                   min_num_of_covering_plies=int(config['CONSTR']['minimum number of covering plies']),
                                   max_num_of_contiguous_plies=int(config['CONSTR']['maximum number of contiguous plies']))


def mat_from_config(config: cp.ConfigParser) -> clt.TransverselyIsotropicMaterial:
    return clt.TransverselyIsotropicMaterial(E1=float(config['PREPREG']['E1']),
                                             E2=float(config['PREPREG']['E2']),
                                             G12=float(config['PREPREG']['G12']),
                                             nu12=float(config['PREPREG']['nu12']),
                                             rho=float(config['PREPREG']['rho']))


def lc_from_config(config: cp.ConfigParser) -> clt.MembraneLoadCase:
    return clt.MembraneLoadCase(Nx=float(config['LOADS']['Nx']),
                                Ny=float(config['LOADS']['Ny']),
                                Nxy=float(config['LOADS']['Nxy']))


def dims_from_config(config: cp.ConfigParser) -> clt.Dimensions:
    return clt.Dimensions(float(config['DIMS']['a']), float(config['DIMS']['b']))


def generic_panel_from_config(config: cp.ConfigParser) -> gs_opt.GenericPanel:
    return gs_opt.GenericPanel(constr=constr_from_config(config),
                               mat=mat_from_config(config),
                               dims=dims_from_config(config))


# For CSV files (given as pandas.DataFrame instances):

# TODO: Implement CSV readers.

def constr_from_df(df: pd.DataFrame) -> list[clt.LaminateConstraints]:
    pass


def mat_from_df(df: pd.DataFrame) -> list[clt.TransverselyIsotropicMaterial]:
    pass


def lc_from_df(df: pd.DataFrame) -> list[clt.MembraneLoadCase]:
    pass


def dims_from_df(df: pd.DataFrame) -> list[clt.Dimensions]:
    pass


# For testing:

if __name__ == '__main__':
    pass
