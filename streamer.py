#!/usr/bin/env python3

# Copyright (c) 2023 Swiss Seismological Service (SED)
# Written by Luca Scarabello @ ETH Zuerich

import os
import sys
import time
import logging
import math
import re
import yaml
from pathlib import Path
from datetime import timedelta
from datetime import datetime
from datetime import timezone
from collections import namedtuple
from obspy.core import UTCDateTime
import obspy as ob

import raw_server as rs


def setup_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%d.%m.%Y %H:%M:%S")
    #
    # create console handler with info log level
    #
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)
    #
    # create file handler which logs even debug messages
    #
    file = logging.FileHandler("streamer.log")
    file.setLevel(logging.INFO)
    file.setFormatter(formatter)
    logger.addHandler(file)
    return logger


#
# Set-up our logger
#
logger = setup_logger("streamer")


def read_file(filename):
    stream = None
    try:
        stream = ob.read(filename)
    except Exception as e:
        logger.error(f"Cannot read {filename}: {e}")
    return stream


class DataChunk:
    def __init__(self, samples, sampling_freq, samples_starttime):
        if samples.ndim != 1:
            raise ValueError("samples is not a 1D array")
        self.samples = samples
        self.sampling_freq = sampling_freq
        self.samples_starttime = samples_starttime

    def start_time(self):
        return self.samples_starttime

    def end_time(self):
        if self.sample_count() == 0:
            raise ValueError("Empty")
        microseconds = 1000000. / self.sampling_freq * \
            (self.sample_count() - 1)
        return self.start_time() + timedelta(microseconds=int(round(microseconds)))

    def empty(self):
        return self.samples.size == 0

    def sample_count(self):
        return self.samples.size

    def data(self):
        return self.samples

    def overlaps(self, starttime, endtime):
        if self.sample_count() == 0:
            raise ValueError("Empty")
        return (self.start_time() <= starttime <= self.end_time()) or \
               (self.start_time() <= endtime <= self.end_time()) or \
               (starttime < self.start_time() and self.end_time() < endtime)

    def sample_index(self, time):
        return ((time - self.start_time()) / timedelta(microseconds=1)
                ) * self.sampling_freq / 1000000.

    def sample_time(self, index):
        microseconds = int(round(index * 1000000. / self.sampling_freq))
        return self.start_time() + timedelta(microseconds=microseconds)

    def extract(self, starttime, endtime):

        if not self.overlaps(starttime, endtime):
            return None

        startoff = math.floor(self.sample_index(starttime))
        endoff = math.ceil(self.sample_index(endtime)) + 1

        if startoff < 0:
            startoff = 0
        if endoff > self.sample_count():
            endoff = self.sample_count()

        # requested data
        slice = DataChunk(self.samples[startoff:endoff],
                          self.sampling_freq,
                          self.sample_time(startoff))

        # drop the requested data from internal buffer
        # pluse anything before starttime
        self.samples = self.samples[endoff:]
        self.samples_starttime = self.sample_time(endoff)

        return slice


class ChannelDataBuffer:
    def __init__(self, data_type, endianness, sampling_freq):
        self.data_type = data_type
        self.sampling_freq = sampling_freq
        self.endianness = endianness
        self.chunks = []

    def empty(self):
        return not self.chunks

    def start_time(self):
        return min([c.start_time() for c in self.chunks], default=None)

    def end_time(self):
        return max([c.end_time() for c in self.chunks], default=None)

    def push(self, trace):
        if trace.data.size == 0:
            return
        if trace.stats.sampling_rate != self.sampling_freq:
            logger.warning(
                f"Drop data {trace.stats.starttime}: sampling freq ({trace.stats.sampling_rate}) differs from configuration ({self.sampling_freq})")
            return
        tr_start = trace.stats.starttime.datetime  # UTCDateTime -> datetime
        tr_end = trace.stats.endtime.datetime     # UTCDateTime -> datetime
        for chunk in self.chunks:
            if chunk.overlaps(tr_start, tr_end):
                logger.warning(
                    f"Data {tr_start}~{tr_end} overlaps with existing data")
        try:
            self.chunks.append(
                DataChunk(trace.data, trace.stats.sampling_rate, tr_start))
        except Exception as err:
            logger.warning(f"Drop data: {err}")

    def pop(self, start, end):
        slice_list = []
        for chunk in self.chunks[:]:
            slice = chunk.extract(start, end)
            if slice is None:
                continue
            slice_list.append(slice)
            if chunk.empty():
                self.chunks.remove(chunk)
        slice_list.sort(key=lambda s: s.start_time())
        return slice_list

    def drop_older_than(self, date):
        for chunk in self.chunks[:]:
            if chunk.end_time() < date:
                self.chunks.remove(chunk)


class DataStreamer:

    def __init__(self,
                 packet_len,
                 simulate_current_time_data,
                 catch_up_mode,
                 buffering_time,
                 channel_settings,
                 streaming_servers):
        self.packet_len = packet_len
        self.simulate_current_time_data = simulate_current_time_data
        self.catch_up_mode = catch_up_mode
        self.buffering_time = buffering_time
        self.streaming_start_time = None
        self.streaming_end_time = None
        self.last_data_streamed_time = None
        self.last_report_time = None
        self.utctime_at_streaming_start = None
        self.ns_at_streaming_start = None
        self.to_catch_up = timedelta(0)
        self.servers = []
        self.buffers = {}
        #
        # Expand freq and format for easy look-up by chan id
        #
        sampling_freq = {}
        data_format = {}
        for cfg in channel_settings:
            freq = cfg["sampling_frequency"]
            frmt = cfg["sample_format"]
            for ch_id in cfg["ids"]:
                if str(ch_id).isdigit():
                    ch_id = int(ch_id)
                    sampling_freq[ch_id] = freq
                    data_format[ch_id] = frmt
                else:
                    min, max = ch_id.split('-')
                    for ch_id in range(int(min), int(max) + 1):
                        sampling_freq[ch_id] = freq
                        data_format[ch_id] = frmt
        #
        # Initialize the streaming servers
        #
        streamers = []
        for server in streaming_servers:
            channels = []
            for ch_id in server['channels']:
                if str(ch_id).isdigit():
                    ch_id = int(ch_id)
                    new_ch = rs.Channel(
                        ch_id, sampling_freq[ch_id], sys.byteorder, data_format[ch_id])
                    channels.append(new_ch)
                else:
                    min, max = ch_id.split('-')
                    for ch_id in range(int(min), int(max) + 1):
                        new_ch = rs.Channel(
                            ch_id, sampling_freq[ch_id], sys.byteorder, data_format[ch_id])
                        channels.append(new_ch)

            streamer = rs.Streamer(
                channels, host=server['host'], port=server['port'])
            streamers.append(streamer)
        self.servers = streamers
        #
        # Initialize channel buffers
        #
        for server in self.servers:
            for ch_id in server.channels:
                self.buffers[ch_id] = ChannelDataBuffer(
                    data_format[ch_id], sys.byteorder, sampling_freq[ch_id])

    def start(self):
        for server in self.servers:
            server.start()

    def stop(self):
        for server in self.servers:
            server.stop()

    def buffer_time_span(self):
        start = self.buffer_start_time()
        if start is None:
            return timedelta(0)
        end = self.buffer_end_time()
        if end is None:
            return timedelta(0)
        return end - start

    def buffer_start_time(self):
        return min([b.start_time()
                   for b in self.buffers.values() if not b.empty()], default=None)

    def buffer_end_time(self):
        return max([b.end_time()
                   for b in self.buffers.values() if not b.empty()], default=None)

    def feed(self, stream):
        #
        # Initialize streaming_start_time and streaming_end_time with the
        # first stream received
        #
        if self.streaming_start_time is None:
            for trace in stream:
                if self.streaming_start_time is None:
                    self.streaming_start_time = trace.stats.starttime.datetime
                elif trace.stats.starttime.datetime < self.streaming_start_time:
                    self.streaming_start_time = trace.stats.starttime.datetime
            self.streaming_end_time = self.streaming_start_time
            logger.info(
                f"DataStreamer initialized: first data time={self.streaming_start_time}")
        #
        # Store the stream into the channel buffers
        #
        for ch_id, trace in enumerate(stream, start=1):
            if trace.stats.endtime.datetime < self.streaming_end_time:
                logger.info(
                    f"Discard data older than last streamed data time (channel {ch_id} trace time {trace.stats.starttime})")
                continue
            for server in self.servers:
                if ch_id in server.channels:
                    logger.debug(
                        f"Feeding channel {ch_id} {trace.stats.starttime.datetime} ~ {trace.stats.endtime.datetime}")
                    self.buffers[ch_id].push(trace)

    def stream(self, max_to_stream=None):
        #
        # If no data has ever been fed just exit
        #
        if self.streaming_start_time is None or self.streaming_end_time is None:
            return
        #
        # Initialization, just once
        #
        if self.last_report_time is None:
            self.last_report_time = self.streaming_start_time

        if self.ns_at_streaming_start is None:
            self.ns_at_streaming_start = time.monotonic_ns()  # int
            self.utctime_at_streaming_start = datetime.now(timezone.utc)
            if self.catch_up_mode:
                self.to_catch_up = self.utctime_at_streaming_start - \
                    self.streaming_start_time.replace(tzinfo=timezone.utc)
                self.to_catch_up -= self.buffering_time
            logger.info(f"DataStreamer: streaming started")
            return

        #
        # Keep track of how much time has elapsed since last stream() call
        #
        elapsed_ns = time.monotonic_ns() - self.ns_at_streaming_start
        elapsed = timedelta(microseconds=math.floor(elapsed_ns / 1000))
        current_time = self.streaming_start_time + elapsed + self.to_catch_up

        #
        # Compute how much time of data we may stream
        #
        if max_to_stream is None:
            max_allowed_end_time = current_time
        else:
            max_allowed_end_time = min(current_time, self.streaming_end_time + max_to_stream)

        #
        # Pass to the streaming server as many data samples as they fit in the elapsed time
        #
        while self.streaming_end_time < max_allowed_end_time:

            stream_start = self.streaming_end_time
            stream_end = stream_start + self.packet_len

            # we streamed all the time we could, so exit now
            if stream_end > max_allowed_end_time:
                break

            #
            # Fetch a packet of data from the buffers and stream it
            #
            for ch_id, buffer in self.buffers.items():
                chunk_list = buffer.pop(stream_start, stream_end)
                if chunk_list:
                    self.last_data_streamed_time = stream_end
                for chunk in chunk_list:
                    #
                    # Finally stream the data to all the servers that accept this channel
                    #
                    for server in self.servers:
                        if ch_id in server.channels:
                            if self.simulate_current_time_data:
                                c_start_time = chunk.start_time() - self.streaming_start_time + \
                                    self.utctime_at_streaming_start
                            else:
                                c_start_time = chunk.start_time()
                            server.feed_data(ch_id, c_start_time, 100, chunk.data())

            self.streaming_end_time = stream_end
        #
        # This should not be necessary, but better be safe
        #
        for buffer in self.buffers.values():
            buffer.drop_older_than(self.streaming_end_time)
        #
        # Periodic status report
        #
        if current_time - self.last_report_time > timedelta(seconds=30):
            buffered_time = self.buffer_time_span()
            logger.info(f"DataStreamer report: elapsed time {elapsed} "
                        f"streaming time {max_allowed_end_time} "
                        f"last streamed data {self.last_data_streamed_time} "
                        f"(delay {(current_time-self.last_data_streamed_time).total_seconds()} sec, "
                        f" buffer {buffered_time.total_seconds()} sec)")
            self.last_report_time = current_time


class PathFilter:
    def __init__(self, pattern):
        self.pattern = pattern
        self.re = re.compile(pattern)

    def filter(self, files):
        valid = [f for f in files if self.re.fullmatch(f)]
        return valid

    def listdir(self, path):
        try:
            all_files = set(os.listdir(path=path))
        except Exception as e:
            logger.error(f"Cannot list {path}: {e}")
            return (set(), set())
        else:
            valid = set(self.filter(all_files))
            invalid = all_files - valid
            return (valid, invalid)


class FileScanner:

    def __init__(self, scan_dir, file_filter,
                 ignore_existing_files, starting_file=None):
        self.scan_dir = scan_dir  # Path
        self.file_filter = file_filter  # PathFilter
        self.backlog = set()
        self.processed = set()
        self.ignored = set()

        if ignore_existing_files:
            self.processed, self.ignored = self.file_filter.listdir(
                self.scan_dir)
            if self.processed:
                logger.debug(f"Ignoring files: {self.processed}")
            if self.ignored:
                logger.debug(f"Unknown files, ignoring: {self.ignored}")
        elif starting_file is not None:
            existing_files, self.ignored = self.file_filter.listdir(
                self.scan_dir)
            if self.ignored:
                logger.debug(f"Unknown files, ignoring: {self.ignored}")
            existing_files = sorted(existing_files)
            while existing_files:
                next_file = existing_files.pop(0)
                if next_file < starting_file:  # previous files
                    logger.debug(f"Ignoring file {next_file}")
                    self.processed.add(next_file)
                elif next_file == starting_file:
                    logger.info(f"Start processing with file {next_file}")
                    self.backlog.add(next_file)
                    break

    def __scan(self):
        #
        # Look for new files in the scan_dir
        #
        valid_files, invalid_files = self.file_filter.listdir(self.scan_dir)
        new_files = valid_files - self.processed

        if new_files:
            logger.info(f"Found new files: {new_files}")
            self.backlog |= new_files

        if invalid_files - self.ignored:
            logger.info(
                f"Unknown files, ignoring: {invalid_files - self.ignored}")
            self.ignored |= invalid_files

    def next_file(self):
        #
        # fetch the next file from the backlog; if empty scan for new files
        #
        if not self.backlog:
            self.__scan()

        if not self.backlog:
            return None

        next_file = sorted(self.backlog)[0]
        self.backlog.remove(next_file)
        self.processed.add(next_file)
        return self.scan_dir / next_file


class DirectoryScanner:
    def __init__(self, scan_dir, subdir_pattern, file_pattern, starting_file):
        self.scan_dir = Path(scan_dir)
        self.subdir_filter = PathFilter(subdir_pattern)
        self.file_filter = PathFilter(file_pattern)
        self.starting_file = starting_file

        self.processed = set()
        self.ignored = set()
        self.file_scanner = None

        #
        # Scan everything, do not skip existing files
        #
        if self.starting_file == "oldest":
            logger.info(f"Looking for the oldest file and start from there")
            self.starting_file = ""
        #
        # Skip existing files and start scanning for new ones
        #
        elif self.starting_file == "latest":
            logger.info(f"Looking for the latest file and start from there")
            self.processed, self.ignored = self.subdir_filter.listdir(
                self.scan_dir)
            if self.ignored:
                logger.info(f"Unknown directories, ignoring: {self.ignored}")
            if self.processed:
                next_dir = sorted(self.processed)[-1]
                to_skip = self.processed - {next_dir}
                if to_skip:
                    logger.debug(f"Ignoring directories: {to_skip}")
                next_dir = self.scan_dir / next_dir
                logger.info(f"Monitoring most recent directory {next_dir}")
                self.file_scanner = FileScanner(
                    next_dir, self.file_filter, True)

        #
        # Start scanning from the specified file onwards
        #
        else:
            logger.info(
                f"Start from the selected file: {self.starting_file}")
            starting_file = Path(self.starting_file).resolve()
            if not starting_file.exists():
                raise ValueError(f"{starting_file} does not exists")
            existing_dirs, self.ignored = self.subdir_filter.listdir(
                self.scan_dir)
            if self.ignored:
                logger.info(f"Unknown directories, ignoring: {self.ignored}")
            existing_dirs = sorted(existing_dirs)
            while existing_dirs:
                next_dir = existing_dirs.pop(0)
                next_dir = self.scan_dir / next_dir
                if next_dir < starting_file.parent:  # previous dates
                    logger.debug(f"Ignoring old directory {next_dir}")
                    self.processed.add(next_dir.name)
                elif next_dir == starting_file.parent:
                    logger.info(
                        f"Monitoring starting from directory {next_dir}")
                    self.processed.add(next_dir.name)
                    self.file_scanner = FileScanner(
                        next_dir, self.file_filter, False, starting_file.name)
                    break

    def __scan(self):
        #
        # Look for new dirs in the scan_dir
        #
        valid_dirs, invalid_dirs = self.subdir_filter.listdir(self.scan_dir)
        new_dirs = valid_dirs - self.processed

        if invalid_dirs - self.ignored:
            logger.info(
                f"Unknown directories, ignoring: {invalid_dirs - self.ignored}")
            self.ignored |= invalid_dirs

        if new_dirs:
            logger.debug(f"Found new directories: {new_dirs}")
            #
            # Load a new FileScanner
            #
            new_dirs = sorted(new_dirs)
            while new_dirs:
                next_dir = new_dirs.pop(0)
                next_dir = self.scan_dir / next_dir
                if self.file_scanner is None or next_dir > self.file_scanner.scan_dir:
                    logger.info(f"Switching to new directory {next_dir}")
                    self.processed.add(next_dir.name)
                    self.file_scanner = FileScanner(
                        next_dir, self.file_filter, False)
                    break
                elif next_dir < self.file_scanner.scan_dir:  # previous dates
                    logger.debug(f"Ignoring old directory {next_dir}")
                    self.processed.add(next_dir.name)

    def next_file(self):
        #
        # Fetch the next file in the current subdir
        #
        if self.file_scanner is not None:
            next_file = self.file_scanner.next_file()
            if next_file is not None:
                return next_file
        #
        # If there are no new files check if there is
        # a new directory that will replace the
        # current self.file_scanner
        #
        self.__scan()
        return None


def run(config, starting_file, catch_up_mode):
    """
    Main acquisition loop, run's until ctrl + c.
    """
    logger.info("Starting with the following settings:")
    logger.info(yaml.dump(config))
    #
    # Set the preferred log levels
    #
    logging.getLogger("raw_api").setLevel(logging.INFO)
    logging.getLogger("raw_server").setLevel(logging.INFO)

    #
    # Read configuration
    #
    packet_len = timedelta(microseconds=config["packet_size_microsec"])
    buffering_time = timedelta(seconds=config['buffering_sec'])
    allow_delayed_data = config['allow_delayed_data']

    #
    # start the data streaming servers
    #
    server = DataStreamer(
        packet_len=packet_len,
        simulate_current_time_data=config["simulate_current_time_data"],
        catch_up_mode=catch_up_mode,
        buffering_time=buffering_time,
        channel_settings=config["channel_settings"],
        streaming_servers=config['streaming_servers']
    )
    server.start()

    #
    # start the directory scanner
    #
    logger.info("Starting data acquisition...")
    dir_scanner = DirectoryScanner(
        scan_dir=config['scan_dir'],
        subdir_pattern=config['subdir_pattern'],
        file_pattern=config['file_pattern'],
        starting_file=starting_file
    )

    try:
        #
        # Run until ctrl + c
        #
        can_stream = False
        while True:

            buffered_time = server.buffer_time_span()
            enough_data_to_stream = buffered_time >= packet_len

            #
            # If the streaming buffer is low on data, feed the next available
            # file to the streaming server
            #
            if not enough_data_to_stream or buffered_time < buffering_time:
                next_file = dir_scanner.next_file()
                if next_file is not None:
                    logger.info(
                        f"Loading file {next_file} (buffer {buffered_time.total_seconds()} sec)")
                    stream = read_file(next_file)
                    if stream is not None:
                        server.feed(stream)

            #
            # if allow_delayed_data is true, stop the streaming when the buffer
            # is low on data, to avoid DataStreamer to create data gaps.
            # We will restart the streaming later, when new files have beed fed to the
            # DataStreamer
            #
            if allow_delayed_data and not enough_data_to_stream:
                if can_stream:
                    logger.info(
                        "No data in buffers: stop streaming and wait for new data.")
                    can_stream = False

            #
            # Stream data
            #
            if can_stream:
                #
                # Compute how much time of data we may stream
                #
                if enough_data_to_stream:
                    # Limit the amount of data to be stramed to give a change to
                    # load and buffer other files before the next call to stream()
                    # This prevents gaps
                    max_to_stream = min(buffered_time, timedelta(seconds=1))
                else:
                    # if there is not enough data to be streamed than the files
                    # have not been generated in time and allow_delayed_data is
                    # disabled.So there is no point in limiting the stream, and
                    # it will generate a gap
                    max_to_stream = None
                #
                # Finally steam
                #
                server.stream(max_to_stream)
            #
            # check if we have buffered enough data, in which case start the streaming
            #
            else:
                if buffered_time >= buffering_time:
                    logger.info(
                        f"Buffering completed: {buffered_time.total_seconds()} sec, start streaming...")
                    can_stream = True

            # we don't want to use 100% cpu, so sleep a while
            sleep_sec = 0.250
            time.sleep(sleep_sec)

    except KeyboardInterrupt:
        logger.info("KeyboardInterrupt detected, exiting...")

    # shutdown
    server.stop()

    logger.info("Exiting, good bye!")


if __name__ == '__main__':

    if len(sys.argv) not in (2, 3, 4):
        print(f"""
Usage:
   {sys.argv[0]} config.yaml [starting-file] [--catch-up]

E.g.
   {sys.argv[0]} config.yaml
   {sys.argv[0]} config.yaml latest
   {sys.argv[0]} config.yaml oldest --catch-up
   {sys.argv[0]} config.yaml scan_dir/yyyy_mm_dd/file.seg2 --catch-up

starting-file (optional):
  Specify from which file to start reading and streaming data. It is optional
  and the default value is "latest"

  Special values are:
   "oldest" : start from the oldest file (the first)
   "latest" : skip all existing files and start from the latest (excluded)


--catch-up (optional):
  Instruct the server to stream the data as fast as possible until it reaches
  the current time (minus buffering_sec). Without this option the data is emitted
  at its original rate, that is 10 seconds of data are streamed in 10 seconds

""")
        sys.exit(0)

    cfg_path = sys.argv[1]
    starting_file = sys.argv[2] if len(sys.argv) >= 3 else "latest"
    catch_up = True if (len(sys.argv) >= 4 and sys.argv[3] == "--catch-up") else False

    config = None
    try:
        with open(cfg_path) as f:
            config = yaml.safe_load(f)
    except Exception as err:
        logger.error(f"Error loading configuration file {cfg_path}: {err}")
        exit(-1)

    run(config, starting_file, catch_up)
