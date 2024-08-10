"""
Contains auxiliary functions meant to be used with pandas.
"""

import csv

import pandas as pd


def prepare_csv_for_pandas(path: str, delimiter: str = ',') -> None:
    """
    Adds delimiters to a CSV such that every row has the same number of cells.
    Necessary because pandas doesn't like rows of different lengths for whatever reason.
    :param path: absolute path to the CSV file
    :param delimiter: CSV delimiter (default: comma)
    :return: nothing
    """
    with open(path, 'r') as file:
        reader = csv.reader(file, delimiter=delimiter)
        rows = list(reader)
    max_delimiters = max(len(row) - 1 for row in rows)
    for row in rows:
        missing_delimiters = max_delimiters - (len(row) - 1)
        row.extend([''] * missing_delimiters)
    with open(path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=delimiter)
        writer.writerows(rows)


def read_csv(path: str, delimiter: str = ',') -> pd.DataFrame:
    prepare_csv_for_pandas(path=path, delimiter=delimiter)
    # noinspection PyTypeChecker
    return pd.read_csv(path, delimiter=delimiter, na_values=None, keep_default_na=False)
