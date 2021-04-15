# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:51:37 2016

@author: ckirst
"""

import numpy as np
import skeleton_graph as sg

g = sg.gt.load_graph('vasculature_graph.gt')

sg.plot_gt_graph_3d(g);