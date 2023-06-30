# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:41:10 2018

@author: Hauke
"""
import numpy as np
import os


class DelPhynMat:

    """
    An Python module to edit DELPHIN m6 files.

    Contains methods to read/write base properties or material functions from/to m6 files

    """

    def __init__(self, fileName):
        """
        :param fileName: the m6 filepath
        """
        self.fileName = fileName


    def __readMaterialFile(self):
        if os.path.isfile(self.fileName):
            f = open(self.fileName, errors="surrogateescape")
            self.fileContent = f.readlines()
            f.close()
        else:
            raise Exception('Couldnt open file: {}'.format(self.fileName))



    def __writeMaterialFile(self):
        if os.path.isfile(self.fileName):
            f = open(self.fileName, 'w')
            f.writelines(self.fileContent)
            f.close()
        else:
            raise Exception('Couldnt open file: {}'.format(self.fileName))


    def __lineOfBasePropertyKeyword(self, keyword):
        keywordList = [keyword,]
        if keyword == 'CE':
            keywordList.append('CET')
        if keyword == 'THETA_POR':
            keywordList.append('OPOR')
        if keyword == 'THETA_CAP':
            keywordList.append('OCAP')
        if keyword == 'THETA_EFF':
            keywordList.append('OEFF')

        for n in range(len(self.fileContent)):
            for key in keywordList:
                if key + '   ' in self.fileContent[n]:
                    return n

        raise Exception("Base property {0:} not found in file '{1:}'".format(keyword, self.fileName))


    def baseProperty(self, keyword, to_float=True):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        ind1 = self.fileContent[line].find('=')+2
        ind2 = self.fileContent[line].find(' ', ind1)
        if to_float:
            return float(self.fileContent[line][ind1:ind2])
        else:
            return self.fileContent[line][ind1:ind2]


    def basePropertyString(self, keyword):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        ind1 = self.fileContent[line].find('=') + 2
        ind2 = self.fileContent[line].find('\n', ind1)
        return self.fileContent[line][ind1:ind2]


    def writeBaseProperty(self, keyword, value):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        ind1 = self.fileContent[line].find('=') + 2
        ind2 = self.fileContent[line].find(' ', ind1)
        self.fileContent[line] = self.fileContent[line][:ind1] + str(value) + self.fileContent[line][ind2:]
        self.__writeMaterialFile()


    def writeBasePropertyString(self, keyword, string):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        ind1 = self.fileContent[line].find('=') + 2
        ind2 = self.fileContent[line].find('\n', ind1)
        self.fileContent[line] = self.fileContent[line][:ind1] + string + self.fileContent[line][ind2:]
        self.__writeMaterialFile()


    def identificationProperty(self, keyword):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        return self.fileContent[line].split('=')[1]


    def writeIdentificationProperty(self, keyword, value):
        self.__readMaterialFile()
        line = self.__lineOfBasePropertyKeyword(keyword)
        ind1 = self.fileContent[line].find('=') + 2
        self.fileContent[line] = self.fileContent[line][:ind1] + str(value)
        self.__writeMaterialFile()


    def materialFunction(self, keyword):
        self.__readMaterialFile()
        for n in range(len(self.fileContent)):
            if keyword in self.fileContent[n] or keyword.replace('Theta_l','Ol') in self.fileContent[n]:
                read = n + 1
                break
        else:
            raise Exception('keyword not found')

        xline = self.fileContent[read].split()
        yline = self.fileContent[read + 1].split()

        xdata = np.zeros(len(xline))
        ydata = np.zeros(len(xline))
        for n in range(len(xline)):
            xdata[n] = float(xline[n])
            ydata[n] = float(yline[n])

        return xdata, ydata


    def writeMaterialFunction(self, xdata, ydata, keyword):
        def floatToLine(x):
            s = '\t'
            for n in x:
                s += '{:9.6}\t'.format(n)
            return s

        if len(xdata) != len(ydata):
            raise Exception('data arrays must be of same length!')

        self.__readMaterialFile()

        # find line to change
        for n in range(len(self.fileContent)):
            if keyword in self.fileContent[n] or keyword.replace('Theta_l','Ol') in self.fileContent[n]:
                change = n + 1
                break
        else:
            raise Exception('keyword {} not found'.format(keyword))
        # change respective lines
        self.fileContent[change] = floatToLine(xdata) + '\n'
        self.fileContent[change + 1] = floatToLine(ydata) + '\n'

        self.__writeMaterialFile()



