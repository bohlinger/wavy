locations = {#'ekofisk':     [56.54, 3.21],
            'ekofiskL':    [56.54, 3.21],
            'gullfaksc':    [61.20, 2.27],
            'asgardb':      [65.11 ,6.79],
            'draugen':      [64.35 ,7.78],
            'gjoa':         [61.33 , 3.9],
            'grane':        [59.17, 2.48], # 2 2 1 
            'heimdal':      [59.57, 2.23], # 2 2 1 
#            'huldra':       [60.87, 2.65], # 2 2 - 
            'kristin':      [65.00, 6.55], # 2 2 2 
            'kvitebjorn':   [61.08, 2.49], # 2 2 1 
#            'njorda':       [64.27, 7.20], # 3 3 1 
            'norne':        [66.02, 8.08], # 3 3 1 
#            'ormenlange':   [63.38, 5.31], # 2 2 1 
            'oseberg':      [60.50, 2.80], # - 2 1 
            'osebergc':     [60.64, 2.78], # 2 2 1 
#            'osebergsyd':   [60.40, 2.80], # 3 2 1 
            'sleipner':     [58.37, 1.91], # 4 4 2 
#            'sleipnerb':    [58.43, 1.70], # 2 - - 
            'snorrea':      [61.45, 2.15], # 4 4 1 
            'snorreb':      [61.52, 2.21], # 2 2 1 
            'stafjorda':    [61.25, 1.85], # - 1 1 
            'statfjordb':   [61.20, 1.83], # 3 3 1 
            'trolla':       [60.64, 3.72], # 3 3 1 
            'trollb':       [60.77, 3.50], # 3 3 - 
            'trollc':       [60.89, 3.60], # 3 1 
            'ula':          [57.10 ,2.85], # 3 3 1 
            'valhall':      [56.28, 3.39], # 5 5 2 
#            'veslefrikka': [60.78, 2.89], # 1 1 1 
            'veslefrikkb':  [60.78, 2.90], # 2 2 2 
            'visund':       [61.37, 2.46], # 4 1 
#            'westnav':     [??? ]
            'heidrun':      [65.32, 7.32],
#            'landegode':    [67.56, 14.17],
             'goliat':      [71.29, 22.32]}

ARCMFClocations = {#'ekofisk':     [56.54, 3.21],
            'ekofiskL':    [56.54, 3.21],
            'gullfaksc':    [61.20, 2.27],
#            'asgardb':      [65.11 ,6.79], corrupf file for Mai 2017
            'draugen':      [64.35 ,7.78],
#            'gjoa':         [61.33 , 3.9], # no file for April 2017
            'grane':        [59.17, 2.48], # 2 2 1 
            'heimdal':      [59.57, 2.23], # 2 2 1 
#            'huldra':       [60.87, 2.65], # 2 2 - 
#            'kristin':      [65.00, 6.55], # no file for Mai 2017 
            'kvitebjorn':   [61.08, 2.49], # 2 2 1 
#            'njorda':       [64.27, 7.20], # 3 3 1 
            'norne':        [66.02, 8.08], # 3 3 1 
#            'ormenlange':   [63.38, 5.31], # 2 2 1 
            'oseberg':      [60.50, 2.80], # - 2 1 
            'osebergc':     [60.64, 2.78], # 2 2 1 
#            'osebergsyd':   [60.40, 2.80], # 3 2 1 
            'sleipner':     [58.37, 1.91], # 4 4 2 
#            'sleipnerb':    [58.43, 1.70], # 2 - - 
            'snorrea':      [61.45, 2.15], # 4 4 1 
            'snorreb':      [61.52, 2.21], # 2 2 1 
            'stafjorda':    [61.25, 1.85], # - 1 1 
#            'statfjordb':   [61.20, 1.83], # 3 3 1 #no file for Mai 2017
            'trolla':       [60.64, 3.72], # 3 3 1 

#            'trollb':       [60.77, 3.50], # no file for April 2017
#            'trollc':       [60.89, 3.60], # no file for April 2017
#            'ula':          [57.10 ,2.85], # 3 3 1 
#            'valhall':      [56.28, 3.39], # 5 5 2 
#            'veslefrikka': [60.78, 2.89], # 1 1 1 
#            'veslefrikkb':  [60.78, 2.90], # no file for Mai 2017
            'visund':       [61.37, 2.46], # 4 1 
#            'westnav':     [??? ]
#            'heidrun':      [65.32, 7.32], # no file for  Mai 2017
#            'landegode':    [67.56, 14.17],
#             'goliat':      [71.29, 22.32] # corrupt file for Mai 2017
                  }


landegode =  {'landegode':    [67.56, 14.17]}

testlocations = {'ekofiskL':     [56.54, 3.21],
                 'draugen':      [64.35 ,7.78]}

bestlocations = testlocations

newlocations = {'goliat':      [71.29, 22.32]}


WMsensors = {'ekofisk':     ['waverider', 'laser altimeter', 'Mir. RangeFinder altimeter'],
             'ekofiskL':    ['waverider', 'laser altimeter', 'Saab WaveRadar altimeter'],
            'gullfaksc':    ['Miros MkIII radar'],
            'asgardb':      ['Miros MkIII radar'],
            'draugen':      ['Miros MkIII radar','Miros MkIII radar'],
            'gjoa':         ['Miros MkIII radar'],
            'grane':        ['Miros MkIII radar'], # 2 2 1 
            'heimdal':      ['Miros MkIII radar'], # 2 2 1 
            'huldra':       [], # 2 2 - 
            'kristin':      [], # 2 2 2 
            'kvitebjorn':   [], # 2 2 1 
            'njorda':       ['Miros MkIII radar'], # 3 3 1 
            'norne':        ['Miros MkIII radar'], # 3 3 1 
            'ormenlange':   ['Miros MkIII radar'], # 2 2 1 
            'oseberg':      ['Miros MkIII radar'], # - 2 1 
            'osebergc':     [], # 2 2 1 
            'osebergsyd':   [], # 3 2 1 
            'sleipner':     ['Miros MkIII radar',''], # 4 4 2 
            'sleipnerb':    [], # 2 - - 
            'snorrea':      ['Miros MkIII radar'], # 4 4 1 
            'snorreb':      [], # 2 2 1 
            'stafjorda':    ['Miros MkIII radar'], # - 1 1 
            'statfjordb':   [], # 3 3 1 
            'trolla':       ['Miros MkIII radar'], # 3 3 1 
            'trollb':       ['Miros MkIII radar'], # 3 3 - 
            'trollc':       [],
            'ula':          [],
            'valhall':      [],
            'veslefrikka':  [],
            'veslefrikkb':  [], 
            'visund':       ['Miros MkIII radar'], 
            'westnav':      ['Miros MkIII radar'],
            'heidrun':      ['Miros MkIII radar','Miros AirGap laser altimeter'], #order unknown
            'goliat':       []} 

bestWMsensor = {'ekofiskL':  2, # python indexing!
            'draugen':      0,
            'kristin':      1,
            'norne':        0,
            'sleipner':     0,
            'valhall':      1,
            'veslefrikkb':  1,
            'goliat':       0}



description = {'ekofisk':     [],
            'gullfaksc':    [],
            'asgardb':      [],
            'draugen':      [],
            'gjoa':         [],
            'grane':        [], # 2 2 1 
            'heimdal':      [], # 2 2 1 
            'huldra':       [], # 2 2 - 
            'kristin':      [], # 2 2 2 
            'kvitebjorn':   [], # 2 2 1 
            'njorda':       [], # 3 3 1 
            'norne':        [], # 3 3 1 
            'ormenlange':   [], # 2 2 1 
            'oseberg':      [], # - 2 1 
            'osebergc':     [], # 2 2 1 
            'osebergsyd':   [], # 3 2 1 
            'sleipner':     [], # 4 4 2 
            'sleipnerb':    [], # 2 - - 
            'snorrea':      [], # 4 4 1 
            'snorreb':      [], # 2 2 1 
            'stafjorda':    [], # - 1 1 
            'statfjordb':   [], # 3 3 1 
            'trolla':       [], # 3 3 1 
            'trollb':       [], # 3 3 - 
            'trollc':       [], # 3 1 
            'ula':          [], # 3 3 1 
            'valhall':      [], # 5 5 2 
            'veslefrikka':  [], # 1 1 1 
            'veslefrikkb':  [], # 2 2 2 
            'visund':       [], # 4 1 
            'heidrun':      [],
            'goliat':       []}


