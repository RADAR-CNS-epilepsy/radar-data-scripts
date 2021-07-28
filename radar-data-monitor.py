#!/usr/bin/env python3

import sys, os, glob
import argparse, fileinput
import csv, gzip
from datetime import datetime, timedelta
import math, re
import numpy as np
from pprint import pprint

E4_SAMPLERATE_ACC = 32
E4_SAMPLERATE_EDA = 4
E4_SAMPLERATE_PPG = 64

PHONE_SAMPLERATE = 5

HOUR = timedelta(hours=1)
SECOND = timedelta(seconds=1)

FILEMASK = "%Y%m%d_%H%M"

E4_TOPICS = {
    "android_empatica_e4_acceleration": E4_SAMPLERATE_ACC,
    "android_empatica_e4_electrodermal_activity": E4_SAMPLERATE_EDA,
    "android_empatica_e4_blood_volume_pulse": E4_SAMPLERATE_PPG
}

PHONE_TOPICS = {
    "android_phone_acceleration": PHONE_SAMPLERATE,
    "android_phone_gyroscope": PHONE_SAMPLERATE,
    "android_phone_magnetic_field": PHONE_SAMPLERATE,
    "android_phone_step_count": PHONE_SAMPLERATE,
    "android_phone_light": PHONE_SAMPLERATE,
    "android_phone_battery_level": PHONE_SAMPLERATE
}

CHECK_TOPICS = {**E4_TOPICS, **PHONE_TOPICS}



if __name__=="__main__":
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass
    cmdline = argparse.ArgumentParser(description="RADAR-CNS data monitor", formatter_class=Formatter)
    cmdline.add_argument('datadirs', type=str, nargs='+', help="data directories\n") # input directories
    cmdline.add_argument('--patient-map', metavar="PATH", type=str, help='Path to CSV file mapping RADAR-IDs to MP-IDs.\n')
    cmdline.add_argument('-l', '--lookback', metavar="HOURS", type=int, default=24, help='Number of hours to look back for coverage report.\nSpecify <1 for infinity\n')
    cmdline.add_argument('--only-active', help='Only process subjects marked as active.\n', action="store_true")
    args = cmdline.parse_args()

    print("Processing directories {}".format(args.datadirs))

    for d in [ os.path.abspath(d) for d in args.datadirs]:
        study_id = os.path.basename(d)
        print()
        print("Study ID: {}".format(study_id))
        print("Data coverage report for the last {} hours:".format(args.lookback))

        if not args.patient_map:
            map_path = os.path.join(d, study_id+".csv")
        else:
            map_path = args.patient_map

        if not os.path.isfile(map_path):
            print("No ID mapping file found for [{}], skipping.".format(d))
            continue
            #raise ValueError("ID mapping file not found!")

        patients = []
        with open(map_path) as map_csvfile:
            map_csvreader = csv.DictReader(map_csvfile)
            patients = [ p for p in map_csvreader ]

        coverage = {} # overall coverage store, one array per pat ID filled with mean coverage per topic
        for p in patients:
            if args.only_active and p["status"] != "active": continue
            coverage[p["RADAR-ID"]] = {} # initialize empty
            for t in CHECK_TOPICS:
                coverage[p["RADAR-ID"]][t] = [np.nan] # initialize empty
                expected_lines = 3600 * CHECK_TOPICS[t] # 1h in sec * sample rate
                # get all the csv files in the topic folder
                inpath = os.path.join(d, p["MP-ID"], t, "*.csv*")
                infiles = sorted(glob.glob(inpath))
                if not infiles: continue

                # get the first and last hour timestamps from the filenames
                f_first = os.path.basename(infiles[0]).split('.')[0]
                f_last = os.path.basename(infiles[-1]).split('.')[0]
                dt_first = datetime.strptime(f_first, FILEMASK)
                dt_last = datetime.strptime(f_last, FILEMASK)

                # create a gapless timeline at 1h interval, initialize each with 0
                if args.lookback < 1:
                    num_hours = int((dt_last-dt_first).total_seconds()/3600)
                    expected_hours = [ dt_first + x * HOUR for x in range(num_hours+1) ]
                else:
                    dt_now = datetime.utcnow().replace(minute=0, second=0, microsecond=0)
                    expected_hours = [ dt_now - x * HOUR for x in range(args.lookback) ]
                topic_coverage = {datetime.strftime(h, FILEMASK):0 for h in expected_hours}

                # iterate data files per topic, fill timeline with data coverage from available files
                for f in infiles:
                    fsplit = os.path.basename(f).split('.')
                    num_lines = 0

                    if fsplit[0] not in topic_coverage.keys(): continue
                    if fsplit[-1] == "gz": num_lines = sum(1 for line in gzip.open(f))
                    elif fsplit[-1] == "csv": num_lines = sum(1 for line in open(f))
                    else: continue

                    file_coverage = num_lines/expected_lines
                    topic_coverage[os.path.basename(f).split('.')[0]] = file_coverage if file_coverage < 1 else 1

                # add mean coverage per topic to array of current pat ID
                cov_vals = np.fromiter(topic_coverage.values(), dtype=float)
                coverage[p["RADAR-ID"]][t] = cov_vals

            #pprint(coverage[p["RADAR-ID"]])
            print("[{}] | {} - {} - ({:8s}) --> ACC:{:<4.0%} EDA:{:<4.0%} BVP:{:<4.0%} PHONE:{:<4.0%}".format(
                    datetime.now().replace(microsecond=0).isoformat(),
                    p["MP-ID"], p["RADAR-ID"], p["status"],
                    np.mean(coverage[p["RADAR-ID"]]["android_empatica_e4_acceleration"]),
                    np.mean(coverage[p["RADAR-ID"]]["android_empatica_e4_electrodermal_activity"]),
                    np.mean(coverage[p["RADAR-ID"]]["android_empatica_e4_blood_volume_pulse"]),
                    np.mean(coverage[p["RADAR-ID"]]["android_phone_acceleration"])
                ))
