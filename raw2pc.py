# this script convert the raw data to pulse compressed data
import KoronaScript.Modules as ksm
import KoronaScript as ks
import os
import glob

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
        print(' ')
        name = channels[channel]['channel_names']
        # just pick the first frequency in the file as the main freq
        MainFrequency = channels[channel]['transducer_frequency'][0] // 1000

        comment = 'Processing pc_' + channel + ' consisting of ' + str(name)
        print(comment)

        # Instantiate the class
        ksi = ks.KoronaScript()
        ksi.add(ksm.Comment(LineBreak='false', Label=comment))
        ksi.add(ksm.ChannelRemoval(Channels=channels[channel]['channels'], KeepSpecified='true'))
        ksi.add(ksm.EmptyPingRemoval())
        ksi.add(ksm.NetcdfWriter(Active="true",
                                 DirName='pc_' + str(channel),
                                 MainFrequency=str(MainFrequency),
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
