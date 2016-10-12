def query_TAP(tap_endpoint, adql_query, table_to_upload=None):
    """
    Query a TAP service (designated by its tap_endpoint)
    with a given ADQL query
    
    Query is performed synchronously
    
    Return an AstroPy Table object
    """
    import requests
    from astropy.table import Table
    from astropy.io.votable import parse_single_table
    import os
    import tempfile
    import warnings
    
    r = requests.post(tap_endpoint + '/sync', data={'query': adql_query, 'request': 'doQuery', 'lang': 'adql', 'format': 'votable', 'phase': 'run'})
    
    #with warnings.catch_warnings():
        #warnings.simplefilter("ignore")
        
    tmp_vot = tempfile.NamedTemporaryFile(delete = False)
    with open(tmp_vot.name, 'w') as h:
        for line in r.iter_lines():
            if line:
                h.write(line.decode(r.encoding)+'\n')

    table = parse_single_table(tmp_vot.name).to_table()

    # finally delete temp files
    os.unlink(tmp_vot.name)

    return table

def plot_scatter_density(data, xcol, ycol, xlabel, ylabel, title, xlim=None, ylim=None):
    import matplotlib.pyplot as plt
    import matplotlib
    #%matplotlib inline
    from scipy.stats import gaussian_kde
    import numpy as np

    x = np.reshape(np.array(data[xcol], copy=False).astype('float'), (len(data['B-V'])))
    y = np.reshape(np.array(data[ycol], copy=False).astype('float'), len(data['B-V']))

    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    
    
    
    w = h = 8
    fig, ax = plt.subplots(figsize = (w, h))
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_title(title)
    ax.scatter(x, y, c=z, s=2, edgecolor='', cmap='jet')
    ax.invert_yaxis()
    plt.show()
