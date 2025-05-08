#!/usr/bin/env python3

import click
import os
import yaml
from wavy.wconfig import load_or_default, load_minimal

@click.command(context_settings={"ignore_unknown_options": True})
@click.option('--f', type=str, default=None,
        help='name of config default file')
@click.option('--path', type=str, default=None,
        help='target path')
@click.option('--t', type=str, default='default',
        help='type (default or minimal)')
def main(f, path, t):
    if path is None:
        path = os.getcwd() 
        print(path)
    fstr = os.path.join(path, f)
    print(fstr)
    if t == 'minimal':
        default_cfg = load_minimal(f)
    elif t == 'default':
        default_cfg = load_or_default(f)

    # Write the dictionary to a YAML file
    with open(fstr, 'w') as file:
        yaml.dump(default_cfg, file, default_flow_style=False,
                  sort_keys=False, width=80)


if __name__ == "__main__":
    main()

