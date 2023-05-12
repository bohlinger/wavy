import wavy.wconfig

def test_load_default():
    c = wavy.wconfig.load_or_default('model_cfg.yaml')

    assert isinstance(c, dict)
    print(c)

    assert c['ARCMFC3']['vardef']['Hs'] == 'VHM0'
