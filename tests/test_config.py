import wavy.wconfig

def test_load_default():
    c = wavy.wconfig.load_or_default('config/model_specs.yaml')

    assert isinstance(c, dict)
    print(c)

    assert c['ARCMFC3']['vars']['Hs'] == 'VHM0'
