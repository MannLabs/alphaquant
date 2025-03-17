import pyarrow
import pandas as pd
import logging
LOGGER = logging.getLogger(__name__)


def read_file(file_path, decimal=".", usecols=None, chunksize=None, sep=None, nrows=None):
    file_path = str(file_path)
    if ".parquet" in file_path:
        if nrows is not None:
            LOGGER.warning(f"nrows parameter is set, but not supported for parquet files. Ignoring nrows parameter.")
        return _read_parquet_file(file_path, usecols=usecols, chunksize=chunksize)
    else:
        if sep is None:
            if ".csv" in file_path:
                sep = ","
            elif ".tsv" in file_path:
                sep = "\t"
            else:
                sep = "\t"
                LOGGER.info(
                    f"neither of the file extensions (.tsv, .csv) detected for file {file_path}! Trying with tab separation. In the case that it fails, please provide the correct file extension"
                )
        return pd.read_csv(
            file_path,
            sep=sep,
            decimal=decimal,
            usecols=usecols,
            encoding="latin1",
            chunksize=chunksize,
            nrows=nrows,
        )


def _read_parquet_file(file_path, usecols=None, chunksize=None):
    if chunksize is not None:
        return _read_parquet_file_chunkwise(
            file_path, usecols=usecols, chunksize=chunksize
        )
    return pd.read_parquet(file_path, columns=usecols)


def _read_parquet_file_chunkwise(file_path, usecols=None, chunksize=None):
    parquet_file = pyarrow.parquet.ParquetFile(file_path)
    for batch in parquet_file.iter_batches(columns=usecols, batch_size=chunksize):
        yield batch.to_pandas()
