#!/usr/bin/env python

from astropy.io import fits
from astropy import units as u
from astropy import wcs 
from astropy.coordinates import SkyCoord
from multiprocessing import Process, Value, Array, Manager
from multiprocessing.pool import ThreadPool
import numpy as np
import sys
import os
import collections
import zerorpc
import logging
import math
import json
import threading
import averager


# import struct

filePrefix = "/data3/artemix-data/FITS/"

class HelloRPC(object):
    
    openedFiles = dict()
    data = dict()
    header = dict()
    averager = dict()
    thread = dict()
    
    def degToHMSDMS(self, ra, dec, fileName, sessionID):
        logging.debug( "in HMSDMS: %d %d" % (ra, dec))
        w = wcs.WCS(self.header[fileName])
        pixcrd = np.array([[ra, dec, 0, 0]], np.float_)
        world = w.wcs_pix2world(pixcrd, 0)
        res = SkyCoord(world[0][0], world[0][1], unit="deg").to_string('hmsdms')
        return res

    def rangeToHMS(self, start, end, step, fileName, sessionID):
        start = int(start)
        end = int(end)
        step = int(step)
        logging.debug( "Getting HMS range.")
        logging.debug("%d %d %d" % (start, end, step))
        res = []
        for i in range(start, end, step):
            logging.debug( "i: %d" % i)
            res.append(self.degToHMSDMS(i, 0, fileName, sessionID).split(" ")[0])
        return res
    
    def rangeToDMS(self, start, end, step, fileName, sessionID):
        start = int(start)
        end = int(end)
        step = int(step)
        logging.info( "Getting DMS range from input parameters '%d %d %d'" % (start, end, step))
        res = []
        for i in range(start, end, step):
            logging.debug( "i: %d" % i)
            res.append(self.degToHMSDMS(0, i, fileName, sessionID).split(" ")[1])
        return res
    
    def _startDataLoading(self, fileName):
        fileNamePath = filePrefix + fileName
        self.averager[fileName] = averager.Averager(str(fileNamePath))
   
   
    def setData(self, fileName, sessionID):
        fileNamePath = filePrefix + fileName
        
        
        
        
        if fileName in self.openedFiles :
            logging.info("File already opened: %s" % fileName)
            self.openedFiles[fileName].add(sessionID)
            d = dict(self.header[fileName])
            return json.dumps(d)
        else :
            logging.debug( "Setting data: %s" % fileName)
            self.openedFiles[fileName] = set()
            self.openedFiles[fileName].add(sessionID)
        
        try:
            hdu_list = fits.open(fileNamePath)
        except IOError as e:
            # logging.info("Error while opening file %s: %s" % (fileName, e))
            raise Exception("Error while opening file %s: %s" % (fileNamePath, e))
        except:
            raise Exception("Unknown error occurred.")
        
        logging.info(fileNamePath)
        self.thread[fileName] = threading.Thread(target=self._startDataLoading, args=(fileName,))
        self.thread[fileName].start()
        
        logging.debug( "Thread started")
      
        self.data[fileName] = hdu_list[0].data[0]
        self.header[fileName] = hdu_list[0].header
        
        self.header[fileName].pop("HISTORY", None) # *None* prevents throwing error
        self.header[fileName].pop("COMMENT", None)
        self.header[fileName].pop("", None)
        try:
            logging.debug("Got RESTFRQ %f" % self.header[fileName]["RESTFRQ"])
        except KeyError:
            self.header[fileName]["RESTFRQ"] = self.header[fileName]["RESTFREQ"]
        
        d = dict(self.header[fileName])
        return json.dumps(d)
        
    def getDimensions(self, fileName, sessionID):
        logging.debug( "Getting dimensions.")
        return self.data[fileName].shape
        
    def flatten(self, l):
        for el in l:
            if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                for sub in self.flatten(el):
                    yield sub
            else:
                yield el
                
    def getSlice(self, sliceNumb, step, fileName, sessionID):
        try:
            sliceNumb = int(sliceNumb)
        except TypeError:
            sliceNumb = 0
      
        logging.info( "Getting slice: %d of %s with shape %r" % (sliceNumb, fileName, self.data[fileName].shape))

        numDimensions = len(self.data[fileName].shape)

        if numDimensions == 2:
            res = self.data[fileName][0:self.data[fileName].shape[0]:step, 0:self.data[fileName].shape[1]:step].tolist()
        elif numDimensions == 3:
            res = self.data[fileName][sliceNumb, 0:self.data[fileName].shape[1]:step, 0:self.data[fileName].shape[2]:step].tolist()
        else:
            logging.debug("Can't process data with such a shape : %r" % self.data[fileName].shape);

        return res
    
    def getFreq(self, fileName, sessionID, x=None, y=None, startZ=None, endZ=None):
        # y x ou x y
        logging.debug( "Getting frequency.")
        if x == None:
            x = 0
        if y == None:
            y = 0
        if startZ == None:
            startZ = 0
        if endZ == None:
            endZ = self.data[fileName].shape[0]
          
        res = self.data[fileName][startZ:endZ, y, x].tolist()
        logging.debug("%d %d" % (x, y))
        return res
        
    def getFreqAverage(self, fileName, sessionID, startY=None, endY=None, startX=None, endX=None):
        logging.debug( "Getting frequency average.")
        
        startZ = 0
        endZ = self.data[fileName].shape[0]
        if startX == None:
            startX = 0
        if endX == None:
            endX = self.data[fileName].shape[2]
        if startY == None:
            startY = 0
        if endY == None:
            endY = self.data[fileName].shape[1]
            
        logging.debug("%d %d %d %d" % (startX, endX, startY, endY))
        logging.debug("%s %s %s %s" % (type(startX), type(endX), type(startY), type(endY)))
        
        # numberElements = (endX - startX) * (endY - startY) 
        res = []
        i = startZ
        pi180 = math.pi / 180 / 4.86e-6 
        #logging.info("observer : " + self.header["OBSERVER"])    

        bmaj = self.header[fileName]["BMAJ"] * pi180
        bmin = self.header[fileName]["BMIN"] * pi180
        convert = math.pi * bmaj * bmin
        cdelt =  4 * math.log(2) * math.fabs(self.header[fileName]["CDELT1"] * pi180) * math.fabs(self.header[fileName]["CDELT2"] * pi180)
        
        while i < endZ: 
            res.append(float(np.nansum(self.data[fileName][i, startY:endY, startX:endX]) / convert * cdelt))
            i += 1
            
        return res
        
    def getAverage(self, fileName, sessionID, step, startZ=None, endZ=None, startY=None, endY=None, startX=None, endX=None):
    
        logging.info("Getting average.")
        
        if startX == None:
            startX = 0
        if endX == None:
            endX = self.data[fileName].shape[2]
        if startY == None:
            startY = 0
        if endY == None:
            endY = self.data[fileName].shape[1]
        if startZ == None:
            startZ = 0
        else:
            startZ = int(startZ)
        if endZ == None:
            endZ = self.data[fileName].shape[0]
        else:
            endZ = int(endZ)
            
        if self.data[fileName].shape[0] == 1:
            res = self.data[fileName].tolist()
            logging.debug("Returning")
            return res

        res = [[]] * endY
        
        # i = startY
        cdelt = 3e+5 * math.fabs(self.header[fileName]["CDELT3"] / self.header[fileName]["RESTFRQ"])
        
        numThreads = 16
        self.thread[fileName].join()
        
        logging.debug("startZ= %d, endZ=%d" % (startZ, endZ))
        self.averager[fileName].setVariables(startX, endX, startY, endY, startZ, endZ, cdelt, numThreads)
        res = self.averager[fileName].countAverage()
        res2 = [list(res[i:i + endX]) for i in range(0, len(res), endX)]
        res3 = np.array(res2)

        logging.info("End of getting average.")
        return res3[0:endY:step, 0:endX:step].tolist()
      
# 
#    End of the class
#
#


def main():
    theRootDir = '/obs/nkasradze/MiningPy'
    theLogDir = theRootDir + '/log'
    theLogFile = theLogDir + '/serverMulti.log'
    logging.basicConfig(filename=theLogFile, format='%(asctime)s %(message)s', level=logging.DEBUG)

    s = zerorpc.Server(HelloRPC())
    s.bind("tcp://*:4244")
    s.run()


if __name__ == '__main__': main()
