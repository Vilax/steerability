#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:22:57 2019

# @author: jlvilas
"""

class emobject:

    def __init__(self, name):
        self.rows = None
        self.cols = None
        self.elems = None
        self.sampling = None
        self.id = None

    def getSampling(self):
        return self.sampling

    def setSampling(self, sampling):
        self.sampling = sampling

    def setEmObjectDimensions(self, obj):
        dims = obj.shape()
        lenDims = len(dims)

        if lenDims == 1:
            self.rows = dims[0]
        if lenDims == 2:
            self.rows = dims[0]
            self.cols = dims[1]
        if lenDims == 3:
            self.rows = dims[0]
            self.cols = dims[1]
            self.elems = dims[2]

    def getDimensions(self):
        return self.rows, self.cols, self.elems

    def geId(self):
        return self.id

    def setId(self, id):
        self.id = id
