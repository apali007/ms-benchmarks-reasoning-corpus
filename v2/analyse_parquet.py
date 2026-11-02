import numpy as np
import pandas as pd
import pyarrow.parquet as pq

pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)
np.set_printoptions(threshold=np.inf)

pf = pq.ParquetFile("/disk/raptor-array/mshab008/msbenchmarking/hardness_index_0/filtered_spectra/hardness_index_0.parquet")
print(pf.schema_arrow)
print(pf.metadata.num_rows)

df = pf.read_row_group(0).to_pandas()
print(df.head(10).T)