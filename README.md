# seedlink_file_streamer
Stream file data (e.g. segy, miniseed) to a seedlink server using the [raw plugin](<https://github.com/swiss-seismological-service/seedlink_raw_plugin>) version v1.1.

This program monitors a directory for new files and when they become available they are streamed to a seedink server. A file can be anything that `obspy.read` can accept (e.g. segy, miniseed, etc). `obspy.read` convert the file to a `Stream` containing multiple `Traces`. Each `Trace` is mapped to a raw_server channel. Read the example `conf.yaml` to understand more about the program behaviour.

Run like this:

```
$ python3 streamer.py

Usage:
   streamer.py config.yaml [starting-file]

starting-file (optional):
  specify from which file to start reading
  and streaming data (e.g. scan_dir/yyyy_mm_dd/file.seg2)

  Special values are:
   "oldest" : start from the oldest file (the first)
   "latest" : skip all existing files and start from the latest (excluded)

  starting-file is optional and the default value is "latest"

```


