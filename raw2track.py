# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import os
from netCDF4 import Dataset
import pandas as pd


"""

This example reads the specified test set (e.g. T2023001), applies pulse compression and stores 
the results as an netcdf. the NetCDF file is read and the pulse compressed data are plotted.

"""


def raw2track(paths, trackingParams):
   # Can paths be set as environment variables???
   # paths = {'inputdir' : pathInputDir
   #          'outputdir' : pathOutputDir
   #          'trranges' : pathToTransducerRanges.xml}
   
   # TransducerRanges.xml contains information on the transducers in the data.
   # Example: 
   """
   <?xml version="1.0" encoding="UTF-8"?>
   <corrections type="RANGE">
      <transducer>
         <parameters>
            <parameter name="Frequency">38</parameter>
            <parameter name="BlindZone">3</parameter>
            <parameter name="Range">35</parameter>
         </parameters>
      </transducer>
      <transducer>
        <parameters>
            <parameter name="Frequency">70</parameter>
            <parameter name="BlindZone">3</parameter>
            <parameter name="Range">35</parameter>
         </parameters>
      </transducer>
      <transducer>
         <parameters>
            <parameter name="Frequency">120</parameter>
            <parameter name="BlindZone">3</parameter>
            <parameter name="Range">35</parameter>
         </parameters>
      </transducer>
      <transducer>
         <parameters>
            <parameter name="Frequency">200</parameter>
            <parameter name="BlindZone">3</parameter>
            <parameter name="Range">35</parameter>
         </parameters>
      </transducer>
      <transducer>
         <parameters>
            <parameter name="Frequency">333</parameter>
            <parameter name="BlindZone">3</parameter>
            <parameter name="Range">35</parameter>
         </parameters>
      </transducer>
   </corrections>
"""   
   # trackingParams can be a dictionary with many values for each key-value pair.  
   # Each key-value pair must have same number of values
   # Example:
   """
   trackingParams = {'Active':                     ["true", "true", "true", "true", "true"],
                 'TrackerType':                    ["Peak", "Peak", "Peak", "Peak", "Peak"],
                 'kHz':                            ["38", "70", "120", "200", "333"],
                 'PlatformMotionType':             ["Floating", "Floating", "Floating", "Floating", "Floating"],
                 'MinTS':                          ["-50","-50","-50","-50","-50"],
                 'PulseLengthDeterminationLevel':  ["50","50","50","50","50"],
                 'MinEchoLength':                  ["0","0","0","0","0"],
                 'MaxEchoLength':                  ["1","1","1","1","1"],
                 'MaxGainCompensation':            ["18","18","18","18","18"],
                 'DoPhaseDeviationCheck':          ["false","false","false","false","false"],
                 'MaxPhaseDevSteps':               ["10","10","10","10","10"],
                 'MaxTS':                          ["0","0","0","0","0"],
                 'MaxDepth':                       ["22","22","22","22","22"], #Must be determined per dataset
                 'MaxAlongshipAngle':              ["10","10","10","10","10"],
                 'MaxAthwartshipAngle':            ["10","10","10","10","10"],
                 'InitiationGateFunction':         [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'InitiationMinLength':            ["1","1","1","1","1"],
                 'GateFunction':                   [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'AlphaBetaEstimator':             [{
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1}
                                                   ],
                 'MaxMissingPings':                ["4","4","4","4","4"],
                 'MaxMissingSamples':              ["2","2","2","2","2"],
                 'MaxMissingPingsFraction':        ["0.7","0.7","0.7","0.7","0.7"],
                 'MinTrackLength':                 ["8","8","8","8","8"],
                 'MinSampleToLengthFraction':      ["0.5","0.5","0.5","0.5","0.5"]}  
   """
   # Example End
   
   # Point to the location of the LSSS installation
   lsss = os.environ["LSSS"]
   ksi = ks.KoronaScript(TransducerRanges = paths['trranges'])
   for ii in range(0, len(trackingParams["kHz"])):
       # Reduce trackingparam dict to only contain the ii-th value in each key-value pair
       reducedTrackingParams = {w:m for w, m in zip(list(trackingParams.keys()), list(list(zip(*list(trackingParams.values())))[ii]))}         
       #add tracking module
       ksi.add(ksm.Tracking(Active =           reducedTrackingParams["Active"],
               TrackerType =                   reducedTrackingParams["TrackerType"],
               kHz=                            reducedTrackingParams["kHz"],
               PlatformMotionType=             reducedTrackingParams["PlatformMotionType"],
               MinTS=                          reducedTrackingParams["MinTS"],
               PulseLengthDeterminationLevel=  reducedTrackingParams["PulseLengthDeterminationLevel"],
               MinEchoLength=                  reducedTrackingParams["MinEchoLength"],
               MaxEchoLength=                  reducedTrackingParams["MaxEchoLength"],
               MaxGainCompensation=            reducedTrackingParams["MaxGainCompensation"],
               DoPhaseDeviationCheck=          reducedTrackingParams["DoPhaseDeviationCheck"],
               MaxPhaseDevSteps=               reducedTrackingParams["MaxPhaseDevSteps"],
               MaxTS=                          reducedTrackingParams["MaxTS"],
               MaxDepth=                       reducedTrackingParams["MaxDepth"], #Must be determined per dataset
               MaxAlongshipAngle=              reducedTrackingParams["MaxAlongshipAngle"],
               MaxAthwartshipAngle=            reducedTrackingParams["MaxAthwartshipAngle"],
               InitiationGateFunction=         reducedTrackingParams["InitiationGateFunction"],
               InitiationMinLength=            reducedTrackingParams["InitiationMinLength"],
               GateFunction=                   reducedTrackingParams["GateFunction"],
               AlphaBetaEstimator=             reducedTrackingParams["AlphaBetaEstimator"],
               MaxMissingPings=                reducedTrackingParams["MaxMissingPings"],
               MaxMissingSamples=              reducedTrackingParams["MaxMissingSamples"],
               MaxMissingPingsFraction=        reducedTrackingParams["MaxMissingPingsFraction"],
               MinTrackLength=                 reducedTrackingParams["MinTrackLength"],
               MinSampleToLengthFraction=      reducedTrackingParams["MinSampleToLengthFraction"]))
       
   
   #Run the script:
   ksi.write()
   ksi.run(src=paths["inputdir"], dst=paths['outputdir'])
   
# Read metadata & env variables
df = pd.read_csv('testdata.csv')
crimac = os.getenv('CRIMACSCRATCH')
# Define input parameters
pathTRanges = "TransducerRanges.xml"
trackingParams = {'Active':                     ["true", "true", "true", "true", "true"],
                 'TrackerType':                    ["Peak", "Peak", "Peak", "Peak", "Peak"],
                 'kHz':                            ["38", "70", "120", "200", "333"],
                 'PlatformMotionType':             ["Floating", "Floating", "Floating", "Floating", "Floating"],
                 'MinTS':                          ["-50","-50","-50","-50","-50"],
                 'PulseLengthDeterminationLevel':  ["50","50","50","50","50"],
                 'MinEchoLength':                  ["0","0","0","0","0"],
                 'MaxEchoLength':                  ["1","1","1","1","1"],
                 'MaxGainCompensation':            ["18","18","18","18","18"],
                 'DoPhaseDeviationCheck':          ["false","false","false","false","false"],
                 'MaxPhaseDevSteps':               ["10","10","10","10","10"],
                 'MaxTS':                          ["0","0","0","0","0"],
                 'MaxDepth':                       ["22","22","22","22","22"], #Must be determined per dataset
                 'MaxAlongshipAngle':              ["10","10","10","10","10"],
                 'MaxAthwartshipAngle':            ["10","10","10","10","10"],
                 'InitiationGateFunction':         [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'InitiationMinLength':            ["1","1","1","1","1"],
                 'GateFunction':                   [{
                                                       "Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20},
                                                       {"Alpha": 2.8,
                                                       "Beta": 2.8,
                                                       "Range": 0.1,
                                                       "TS": 20}],
                 'AlphaBetaEstimator':             [{
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1},
                                                       {
                                                       "Alpha": 0.9,
                                                       "Beta": 0.1}
                                                   ],
                 'MaxMissingPings':                ["4","4","4","4","4"],
                 'MaxMissingSamples':              ["2","2","2","2","2"],
                 'MaxMissingPingsFraction':        ["0.7","0.7","0.7","0.7","0.7"],
                 'MinTrackLength':                 ["8","8","8","8","8"],
                 'MinSampleToLengthFraction':      ["0.5","0.5","0.5","0.5","0.5"]}

for _dataset in df['dataset']:
   print(_dataset)
   inputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                           _dataset, 'ACOUSTIC',
                           'EK80', 'EK80_RAWDATA')
   outputdir = os.path.join(crimac, 'CRIMAC-FM-testdata', _dataset[1:5],
                            _dataset, 'ACOUSTIC',
                            'LSSS', 'KORONA')
   
   if os.path.exists(inputdir):
       print(inputdir)
       print(outputdir)
       #raw2pc(inputdir, outputdir)
       paths = {'inputdir' : inputdir,
               'outputdir' : outputdir,
               'trranges' : pathTRanges}
       raw2track(paths, trackingParams)
