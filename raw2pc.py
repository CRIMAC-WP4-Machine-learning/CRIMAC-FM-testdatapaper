# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import os
import glob
import yaml

import sys
import raw2meta

"""
This example reads the specified test set (e.g. T2023001), applies
pulse compression and stores the results as an netcdf. the NetCDF file
is read and the pulse compressed data are plotted.
"""

def raw2pc(inputdir, outputdir, channels, debug=False):
    """
    Raw2pc convert the raw files to pulse compressed files (when applicable)
    for each ping group using korona and KoronaScript.
    """
    # Loop over the different ping groups
    for channel in channels:
        name = channels[channel]['channel_names']
        # This is only needed for GRIDDED
        # MainFrequency = channels[channel]['transducer_frequency'][0] // 1000

        comment = 'Processing pc_' + channel + ' consisting of ' + str(name)
        print(comment)

        # Instantiate the class
        ksi = ks.KoronaScript()
        ksi.add(ksm.Comment(LineBreak='false', Label=comment))
        ksi.add(ksm.ChannelRemoval(Channels=channels[channel]['channels'], KeepSpecified='true'))
        ksi.add(ksm.EmptyPingRemoval())
        ksi.add(ksm.NetcdfWriter(Active="true",
                                 DirName='pc_' + str(channel),
                                 # MainFrequency=str(MainFrequency),
                                 MaxRange=400,
                                 WriterType="CHANNEL_GROUPS",
                                 GriddedOutputType="PULSE_COMPRESSION",
                                 WriteAngels="true",
                                 FftWindowSize="2",
                                 DeltaFrequency="1",
                                 ChannelGroupOutputType="PULSE_COMPRESSION"))

        if debug: ksi.write()
        ksi.run(src=inputdir, dst=outputdir, debug=debug)

        # Remove temporary korona files
        for f in glob.glob(outputdir + '/*korona.*'): os.remove(f)


if __name__ == '__main__':
    if not len(sys.argv) == 3:
        print(f'Usage: {sys.argv[0]} <inputdir> <outputdir>')
        exit(-1)
    indir, outd = sys.argv[1], sys.argv[2]
    if os.path.exists(outd):
        print(f'Output dir "{outd}" already exists. Aborting.')
        exit(-1)

    os.makedirs(outd, exist_ok=True)
    channels, con, ind = raw2meta.raw2meta(indir)
    print(f'Channels:\n{yaml.dump(channels)}')
    raw2pc(indir, outd, channels, debug=False)

