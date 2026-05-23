# Data model for TS(f)

This documents the data model for target strength by frequency. This is an Xarray data structure, stored in `NetCDF` or `zarray` formats.

```
Groups: Environment, frequency_0-N

Dimension:
frequency                        (frequency) float64
i                                (i) int64

Variables:
compensated_TS                   (i, frequency) float32
single_target_alongship_angle    (i) float32
single_target_athwartship_angle  (i) float32
ping_time                        (i) int64
ping_number                      (i) float32
single_target_range              (i) float32
single_target_identifier         (i) int64
single_target_start_range        (i) float32
single_target_stop_range         (i) float32

Attributes:
single_target_detection_algorithm 
exportTime:                      2024-01-10T12:37:54Z
lsssVersion:                     2.15.0
formatLastChangedInLsssVersion:  2.3.0
parameters:                      {'AllFrequencies': False, 'TargetExtentM...
units:                           {'FrequencyResolution': 'Hz', 'Frequency...
comments:                        ['The format of and definitions in this ...
```