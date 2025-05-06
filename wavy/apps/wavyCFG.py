#!/usr/bin/env python3

import click
import os
import yaml
from wavy.wconfig import load_or_default

@click.command(context_settings={"ignore_unknown_options": True})
@click.option('--f', type=str, default=None,
        help='name of config default file')
@click.option('--path', type=str, default=None,
        help='target path')
def main(f, path):
    if path is None:
        path = os.getcwd() 
        print(path)
    fstr = os.path.join(path, f)
    print(fstr)
    default_cfg = load_or_default(f)

    # Write the dictionary to a YAML file
    with open(fstr, 'w') as file:
        yaml.dump(default_cfg, file, default_flow_style=False)


if __name__ == "__main__":
    main()

