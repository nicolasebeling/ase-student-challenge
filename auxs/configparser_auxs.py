"""
Contains auxiliary functions meant to be used with configparser.
"""

import configparser as cp


def read_config(paths: str | list[str]) -> cp.ConfigParser:
    config = cp.ConfigParser()
    read = config.read(paths)
    output = 'Read the following config files: '
    for i, file in enumerate(read):
        output += file
        if i < len(read) - 1:
            output += ', '
    print(f'{output}\n')
    return config
