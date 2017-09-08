# MiningPy

ZeroRPC is a light-weight, reliable and language-agnostic library for distributed communication between server-side processes. It builds on top of ZeroMQ and MessagePack. Support for streamed responses - similar to python generators - makes zerorpc more than a typical RPC engine. Built-in heartbeats and timeouts detect and recover from failed requests. Introspective capabilities, first-class exceptions and the command-line utility make debugging easy.

ZeroRPC server to calculate average of Radio images

serverMulti.py --- ZeroRPC server to calculate average of Radio images for MiningPlot. 

Calculation of average is done with C++, using SWIG to conenct python and C++;
