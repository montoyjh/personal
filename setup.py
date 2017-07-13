#!/usr/bin/env python

import os

from setuptools import setup

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='personal',
        version='0.1',
        description='Personal sandbox for Joseph Montoya',
    )
