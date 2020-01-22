#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 13:46:09 2018

@author: herminio
"""
import json


class Properties(object):
    def __init__(self, fileName):
        self.fileName = fileName
        self.__readJsonFile()
        self.__getRegions()
        self.numberOfRegions = len(self.regions)

    def fromRegionGetProperty(self, region, prop):
        if type(region) == int:
            try:
                region = self.regions[region]
            except:
                raise Exception('Region %i out of bounds.'%region)

        if self.regions.count(region) != 0:
            props = [p for p in self.data.get(region).keys()]
            if props.count(prop) != 0:
                KEY, = self.data.get(region).get(prop).items()
                return self.data.get(region).get(prop).get(KEY[0])
            else:
                raise Exception('Property %s does not belong to region %s'%(prop, region))
        else:
            raise Exception('Region %s is not specified in file %s'%(region, self.fileName))

    def __getRegions(self):
        self.regions = [region for region in self.data.keys()]

    def __readJsonFile(self):
        f = open(self.fileName, "r")
        self.data = json.load(f)







if __name__ == '__main__':
    fileName = "Json_Files/solid.json"
    p = Properties(fileName)
    print p.data, '\n'
    print p.materials
    print p.fromMaterialGetProperty(p.materials[0], 'Density')
    print p.fromMaterialGetProperty(1, 'Density')
