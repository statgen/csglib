import pandas as pd

def thin_dframe(dframe,columns,digits=2,only_cols=True):
  """
  Given a data frame and column names, thin the data frame down by rounding values
  in the columns to N digits, and then finding the unique rows (using only those columns). 
  The original values are returned, not the rounded values. 
  
  :param dframe: pandas data frame
  :param columns: list of column names
  :param digits: number of digits at which to round
  :param only_cols: should we return the dframe subsetted to columns? 
  """
  
  rounded = []
  for col in columns:
    rounded.append(dframe[col].round(digits))
    
  tmpdf = pd.DataFrame({col : rounded[i] for i, col in enumerate(columns)})
  thin = ~ tmpdf.duplicated()
  
  if only_cols:
    fcols = columns
  else:
    fcols = dframe.columns
  
  return dframe.loc[thin,fcols]
