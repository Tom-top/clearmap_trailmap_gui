#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:32:25 2020

@author: thomas.topilko
"""

import xml.etree.ElementTree as ET;
import xml.dom.minidom;
import os;

svgFile = "/home/thomas.topilko/Desktop/Line_Test.svg";
namespace = "{http://www.w3.org/2000/svg}";

# for parsing svg from a file:
svg = ET.parse(svgFile);

root = svg.getroot();

#for child in root :
#    print("\n")
#    print(child.tag, child.attrib)
    
graphicalObj = root.findall("{0}g".format(namespace))[0];

for child in graphicalObj :
    print("\n")
    print(child.tag, child.attrib)

element = "path";

for child in graphicalObj :
    
    print("\n")
    
    if element in child.tag :
        
        print(child.attrib["id"])
        if child.attrib["id"] == "EW-POA" :
        
            newParams = [];
            
            for param in child.attrib["style"].split(";") :
                
                if "stroke-width" in param :
                    
                    temp = "stroke-width:{0}".format(5);
                    newParams.append(temp);
                
                else :
                    
                    newParams.append(param);
                    
            child.attrib["style"] = ";".join(newParams);
        
svg.write(os.path.join("/home/thomas.topilko/Desktop", "output.svg"))
        
            