#!/usr/bin/env python

from astropy.io import fits
from astropy import units as u
from astropy import wcs 
from astropy.coordinates import SkyCoord
from multiprocessing import Process, Value, Array, Manager
from multiprocessing.pool import ThreadPool
import numpy as np
import sys
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
    data = []
    header = []
    averager = None
    thread = None
    
    def degToHMSDMS(self, ra, dec):
        logging.debug( "in HMSDMS: %d %d" % (ra, dec))
        w = wcs.WCS(self.header)
        pixcrd = np.array([[ra, dec, 0, 0]], np.float_)
        world = w.wcs_pix2world(pixcrd, 0)
        res = SkyCoord(world[0][0], world[0][1], unit="deg").to_string('hmsdms')
        return res

    def rangeToHMS(self, start, end, step):
        start = int(start)
        end = int(end)
        step = int(step)
        logging.debug( "Getting HMS range.")
        logging.debug("%d %d %d" % (start, end, step))
        res = []
        for i in range(start, end, step):
            logging.debug( "i: %d" % i)
            res.append(self.degToHMSDMS(i, 0).split(" ")[0])
        return res
    
    def rangeToDMS(self, start, end, step):
        start = int(start)
        end = int(end)
        step = int(step)
        logging.debug( "Getting DMS range.")
        logging.debug("%d %d %d" % (start, end, step))
        res = []
        for i in range(start, end, step):
            logging.debug( "i: %d" % i)
            res.append(self.degToHMSDMS(0, i).split(" ")[1])
        return res
    
    def _startDataLoading(self, fileName):
        self.averager = averager.Averager(str(fileName))
   
   
    def setData(self, fileName):
        fileName = filePrefix + fileName
        
        
        logging.debug( "Setting data: %s" % fileName)
        
        try:
            hdu_list = fits.open(fileName)
        except IOError as e:
            # logging.info("Error while opening file %s: %s" % (fileName, e))
            raise Exception("Error while opening file %s: %s" % (fileName, e))
        except:
            raise Exception("Unknown error occurred.")
        
        self.thread = threading.Thread(target=self._startDataLoading, args=(fileName,))
        self.thread.start()
        
        logging.debug( "Thread started")
      
        self.data = hdu_list[0].data[0]
        self.header = hdu_list[0].header
        
        self.header.pop("HISTORY", None) # *None* prevents throwing error
        self.header.pop("COMMENT", None)
        self.header.pop("", None)
        try:
            logging.debug("Got RESTFRQ %f" % self.header["RESTFRQ"])
        except KeyError:
            self.header["RESTFRQ"] = self.header["RESTFREQ"]
        
        d = dict(self.header)
        return json.dumps(d)
        
    def getDimensions(self):
        logging.debug( "Getting dimensions.")
        return self.data.shape
        
    def flatten(self, l):
        for el in l:
            if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                for sub in self.flatten(el):
                    yield sub
            else:
                yield el
                
    def getSlice(self, sliceNumb, step):
        try:
            sliceNumb = int(sliceNumb)
        except TypeError:
            sliceNumb = 0
      
        logging.debug( "Getting slice: %d" % sliceNumb)
        res = self.data[sliceNumb, 0:self.data.shape[1]:step, 0:self.data.shape[2]:step].tolist()
        return res
    
    def getFreq(self, y=None, x=None, startZ=None, endZ=None):
        logging.debug( "Getting frequency.")
        if x == None:
            x = 0
        if y == None:
            y = 0
        if startZ == None:
            startZ = 0
        if endZ == None:
            endZ = self.data.shape[0]
          
        res = self.data[startZ:endZ, y, x].tolist()
        logging.debug("%d %d" % (x, y))
        return res
        
    def getFreqAverage(self,  startY=None, endY=None, startX=None, endX=None):
        logging.debug( "Getting frequency average.")
        
        startZ = 0
        endZ = self.data.shape[0]
        if startX == None:
            startX = 0
        if endX == None:
            endX = self.data.shape[2]
        if startY == None:
            startY = 0
        if endY == None:
            endY = self.data.shape[1]
            
        logging.debug("%d %d %d %d" % (startX, endX, startY, endY))
        logging.debug("%s %s %s %s" % (type(startX), type(endX), type(startY), type(endY)))
        
        # numberElements = (endX - startX) * (endY - startY) 
        res = []
        i = startZ
        pi180 = math.pi / 180 / 4.86e-6 
        #logging.info("observer : " + self.header["OBSERVER"])	

	bmaj = self.header["BMAJ"] * pi180
        bmin = self.header["BMIN"] * pi180
        convert = math.pi * bmaj * bmin
        cdelt =  4 * math.log(2) * math.fabs(self.header["CDELT1"] * pi180) * math.fabs(self.header["CDELT2"] * pi180)
        
        while i < endZ: 
            res.append(float(np.nansum(self.data[i, startY:endY, startX:endX]) / convert * cdelt))
            i += 1
            
        return res
        
    def getAverage(self, step, startZ=None, endZ=None, startY=None, endY=None, startX=None, endX=None):
    
        logging.info("Getting average.")
    
        if startX == None:
            startX = 0
        if endX == None:
            endX = self.data.shape[2]
        if startY == None:
            startY = 0
        if endY == None:
            endY = self.data.shape[1]
        if startZ == None:
            startZ = 0
        else:
            startZ = int(startZ)
        if endZ == None:
            endZ = self.data.shape[0]
        else:
            endZ = int(endZ)
            
        if self.data.shape[0] == 1:
            res = self.data.tolist()
            logging.debug("Returning")
            return res

        res = [[]] * endY
        
        # i = startY
        cdelt = 3e+5 * math.fabs(self.header["CDELT3"] / self.header["RESTFRQ"])
        
        numThreads = 16
        self.thread.join()
        
        logging.debug("startZ= %d, endZ=%d" % (startZ, endZ))
        self.averager.setVariables(startX, endX, startY, endY, startZ, endZ, cdelt, numThreads)
        res = self.averager.countAverage()
        res2 = [list(res[i:i + endX]) for i in range(0, len(res), endX)]
        res3 = np.array(res2)

        logging.info("End of getting average.")
        return res3[0:endY:step, 0:endX:step].tolist()
      
# 
#    End of the class
#
#


def main():
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.DEBUG)

    s = zerorpc.Server(HelloRPC())
    s.bind("tcp://*:4242")
    s.run()


if __name__ == '__main__': main()
