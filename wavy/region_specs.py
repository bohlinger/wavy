''' 
dict for regions, currently defined by lat(e/w) lon(s/n) borders.
llcrnrlon=12.2,llcrnrlat=67.6,urcrnrlon=13.2,urcrnrlat=67.9
Ultimately, a shape file or polygon would be preferred.
Polygones are given as dict below.
'''
region_dict = {"Global":            {
                                    "llcrnrlon":-180,
                                    "llcrnrlat":-90,
                                    "urcrnrlon":180,
                                    "urcrnrlat":90
                                    },
                "ARCMFC":           {
                                    "boundinglat":50,
                                    },
                "Arctic":           {
                                    "boundinglat":66,
                                    },
                "Mosk_dom":         {
                                    "llcrnrlon":9.4,
                                    "llcrnrlat":66.,
                                    "urcrnrlon":30,
                                    "urcrnrlat":73.
                                    },
                "Moskenes":         {
                                    "llcrnrlon":12.,
                                    "llcrnrlat":67.5,
                                    #"urcrnrlon":13.5,
                                    "urcrnrlon":14,
                                    #"urcrnrlat":68},
                                    "urcrnrlat":68.5
                                    },
                "MoskWC":         {
                                    "llcrnrlon":12.,
                                    "llcrnrlat":67.5,
                                    #"urcrnrlon":13.5,
                                    "urcrnrlon":14,
                                    #"urcrnrlat":68},
                                    "urcrnrlat":68.5
                                    },
                "MoskNC":         {
                                    "llcrnrlon":12.,
                                    "llcrnrlat":67.5,
                                    #"urcrnrlon":13.5,
                                    "urcrnrlon":14,
                                    #"urcrnrlat":68},
                                    "urcrnrlat":68.5
                                    },
                "Sulafj":           {
                                    "llcrnrlon":5.,
                                    "llcrnrlat":62.,
                                    "urcrnrlon":8.14,
                                    "urcrnrlat":63.15
                                    },
                "mwam4":            {
                                    "llcrnrlon":-32.,
                                    "llcrnrlat":42,
                                    "urcrnrlon":94,
                                    "urcrnrlat":84
                                    },
                "mwam8":            { # problem with plotting
                                    "llcrnrlon":-220.,
                                    "llcrnrlat":36,
                                    "urcrnrlon":140,
                                    "urcrnrlat":84
                                    },
                "ecwam":            {
                                    "llcrnrlon":-180,
                                    "llcrnrlat":-90,
                                    "urcrnrlon":180,
                                    "urcrnrlat":90
                                    },
                "man":              { # problem with plotting
                                    "llcrnrlon":-220.,
                                    "llcrnrlat":36,
                                    "urcrnrlon":140,
                                    "urcrnrlat":84
                                    },
                "Erin1W":           {
                                    "llcrnrlon":-20.,
                                    "llcrnrlat":60.,
                                    "urcrnrlon":80.,
                                    "urcrnrlat":89.
                                    },
                "Erin2W":           {
                                    "llcrnrlon":-20.,
                                    "llcrnrlat":60.,
                                    "urcrnrlon":80.,
                                    "urcrnrlat":89.
                                    }
                }

poly_dict = {"NordicSeas": {
                            "lats":[61.,66.,79.,76.,70.],
                            "lons":[5.,-35.,-19.,15.,19.]
                            },
             "BarentsSea": {None
                            }, 
             "mwam800c3":  {"lats":[58.16,63.35,63.11,57.55],
                           "lons":[2.47,4.56,8.55,6.10]
                            },
             "mwam4":      {"lats":[42.55,53.09,83.06,61.56],
                           "lons":[2.18,30.58,92.53,65.33]
                            },
             "mwam8":      {"lats":[42.18,53.10,36.14,53.16,51.02,54.50,56.33],
                           "lons":[3.21,35.15,91.36,134.23,176.43,160.11,101.23]
                            },
             "arcmfc":     {"lats":[42.18,53.10,36.14,53.16,51.02,54.50,56.33],
                           "lons":[3.21,35.15,91.36,134.23,176.43,160.11,101.23]
                            }
            }
#from matplotlib.patches import Polygon
#from matplotlib.path import Path
#import numpy as np
#poly = Polygon(list(zip(poly_dict['NordicSeas']['lons'], poly_dict['NordicSeas']['lats'])), closed=True)
#lon=np.array([1,2]).ravel()
#lat=np.array([69,70]).ravel()
#points=np.c_[lon,lat]
#Path(poly.xy).contains_points(points)
