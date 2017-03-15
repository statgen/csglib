import codecs
from IPython.display import display, HTML
from six.moves import cStringIO as StringIO

def in_notebook():
  """ 
  Determine if the code is currently executing as a jupyter notebook. 
  """

  return 'ipykernel' in sys.modules

def show_table(df,*args,**kwargs):
  return HTML(df.to_html(*args,**kwargs))

def text(string,*args,**kwargs):
  """
  Helper function to create text in the IPython notebook, using string
  formatting operations. 
  
  Displays the text as an HTML element inline. 

  Example: 
    text("There were {neqtls} eQTL variants discovered.",neqtls=100)
  """

  string = string.format(*args,**kwargs)
  display(HTML("<p>%s</p>" % string))

def df_to_uri(dframe,link_text="Download",filename="table.tab",sep="\t",na_rep="NA",encoding="utf-8",*args,**kwargs):
  """
  Converts a pandas data frame into a data URI HTML element. 

  This returns an HTML <a> element inline in the notebook, with `link_text` shown. 
  When the user clicks on the link, the browser will attempt to download the data frame
  named as `filename`. 
  
  This currently does not work for very large data URIs in Chrome only. 
  See: http://stackoverflow.com/questions/16761927/aw-snap-when-data-uri-is-too-large
  """
  
  buf = StringIO()
  dframe.to_csv(buf,sep=sep,na_rep=na_rep,encoding=encoding,*args,**kwargs)
  
  s = buf.getvalue()
  b64 = codecs.encode(s.encode("utf-8"),"base64")
  
  html = '<a href="data:text/csv;base64,{}" download="{}">{}</a>'.format(b64.decode("utf-8"),filename,link_text)
  return HTML(html)

