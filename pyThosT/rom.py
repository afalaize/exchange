# -*- coding: utf-8 -*-
"""
Created on Thu Feb 02 11:53:24 2017

@author: afalaize
"""

from pod import POD

class ReducedOrderModel:
    def __init__(self, ts_velocity, ts_density, ts_viscosity):
        self.pod = POD(ts_velocity)
        