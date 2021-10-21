#!/usr/bin/env python3

import sys, os, glob, copy
import argparse, fileinput
import csv, gzip, json
from datetime import datetime, timedelta
import pytz
import math, re
import numpy as np
from pprint import pprint

E4_SAMPLERATE_ACC = 32
E4_SAMPLERATE_EDA = 4
E4_SAMPLERATE_PPG = 64

BIOV_SAMPLERATE_ACC = 51.2
BIOV_SAMPLERATE_EDA = 1
BIOV_SAMPLERATE_PPG = 51.2

PHONE_SAMPLERATE = 5

HOUR = timedelta(hours=1)
SECOND = timedelta(seconds=1)

FILEMASK = "%Y%m%d_%H%M"

E4_TOPICS = {
    "android_empatica_e4_acceleration": E4_SAMPLERATE_ACC,
    "android_empatica_e4_electrodermal_activity": E4_SAMPLERATE_EDA,
    "android_empatica_e4_blood_volume_pulse": E4_SAMPLERATE_PPG
}

BIOV_TOPICS = {
    "android_biovotion_vsm1_acceleration": BIOV_SAMPLERATE_ACC,
    "android_biovotion_vsm1_galvanic_skin_response": BIOV_SAMPLERATE_EDA,
    "android_biovotion_vsm1_ppg_raw": BIOV_SAMPLERATE_PPG
}

PHONE_TOPICS = {
    "android_phone_acceleration": PHONE_SAMPLERATE,
    # "android_phone_gyroscope": PHONE_SAMPLERATE,
    # "android_phone_magnetic_field": PHONE_SAMPLERATE,
    # "android_phone_step_count": PHONE_SAMPLERATE,
    # "android_phone_light": PHONE_SAMPLERATE,
    # "android_phone_battery_level": PHONE_SAMPLERATE
}

# https://matplotlib.org/stable/api/markers_api.html
# https://matplotlib.org/stable/gallery/color/named_colors.html
QUESTIONNAIRE_TOPICS = {
    "questionnaire_epi_evening_questionnaire": ('o','cyan'),
    "questionnaire_epi_seizure_diary": ('*','yellow'),
}

COVERAGE_TOPICS = {**E4_TOPICS, **PHONE_TOPICS}

STR_REPLACEMENTS = {
    "android_phone": "PH",
    "android_empatica_e4": "E4",
    "android_biovotion_vsm1": "BV",
    "acceleration": "ACC",
    "electrodermal_activity": "EDA",
    "galvanic_skin_response": "EDA",
    "blood_volume_pulse": "PPG",
    "ppg_raw": "PPG"
}

def getShortTopicStr(topic):
    shortTopic = copy.deepcopy(topic)
    for longStr, shortStr in STR_REPLACEMENTS.items():
        shortTopic = shortTopic.replace(longStr, shortStr)
    return shortTopic



def getFirstLastDay(topics, inpath):
    # get first and last days of recordings for each topic
    first_days = []
    last_days = []
    for t in topics:
        # get all the csv files in the topic folder
        infiles = sorted(glob.glob(os.path.join(inpath, t, "*.csv*")))
        if not infiles: continue

        # get the first and last hour timestamps from the filenames
        f_first = os.path.basename(infiles[0]).split('.')[0][:13]
        f_last = os.path.basename(infiles[-1]).split('.')[0][:13]
        dt_first = datetime.strptime(f_first, FILEMASK)
        dt_last = datetime.strptime(f_last, FILEMASK)
        first_days.append(dt_first.replace(hour=0))
        last_days.append(dt_last.replace(hour=23))

    # return the outer bounds of first and last day over all topics
    return (min(first_days, default=None), max(last_days, default=None))


# https://stackoverflow.com/a/9631635
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


def getTopicCoverage(inpath, expected_lines, dt_first_day, dt_last_day):
    # init, return if no files found
    topic_coverage_all = [np.nan]
    topic_num_samples_all = 0
    infiles = sorted(glob.glob(inpath))
    if not infiles: return (topic_coverage_all, topic_num_samples_all)

    # create the coverage value timelines
    dt_today = datetime.utcnow().replace(minute=0, second=0, microsecond=0, hour=23)
    num_hours_all = int((dt_last_day-dt_first_day).total_seconds()/3600)
    expected_hours_all = [ dt_first_day + x * HOUR for x in range(num_hours_all+1) ]
    topic_coverage_all = {datetime.strftime(h, FILEMASK):0 for h in expected_hours_all}

    # iterate data files per topic, fill timeline with data coverage from available files
    for fname in infiles:
        fsplit = os.path.basename(fname).split('.')
        num_lines = 0

        # read lines, skip if not in timeline
        if fsplit[0] not in topic_coverage_all.keys(): continue
        if fsplit[-1] == "gz":
            # num_lines = sum(1 for line in gzip.open(fname))
            with gzip.open(fname, "rt", encoding="utf-8", errors='ignore') as f:
                num_lines = sum(bl.count("\n") for bl in blocks(f))
        elif fsplit[-1] == "csv":
            # num_lines = sum(1 for line in open(fname))
            with open(fname, "rt", encoding="utf-8", errors='ignore') as f:
                num_lines = sum(bl.count("\n") for bl in blocks(f))
        else: continue

        # calculate coverage, fill in timeline, upper bound at 100%
        file_coverage = num_lines/expected_lines
        if fsplit[0] in topic_coverage_all.keys(): topic_coverage_all[fsplit[0]] = file_coverage if file_coverage < 1 else 1
        topic_num_samples_all += num_lines
    
    return (topic_coverage_all, topic_num_samples_all)




if __name__=="__main__":
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass
    cmdline = argparse.ArgumentParser(description="RADAR-CNS data overview", formatter_class=Formatter)
    cmdline.add_argument('datadirs', type=str, nargs='+', help="data directories\n") # input directories
    cmdline.add_argument('--patient-ids', type=str, nargs='*', help='List of patient IDs (i.e. folder names), takes precedent over \'--patient-map\'.\n')
    cmdline.add_argument('--patient-map', metavar="PATH", type=str, help='Path to CSV file mapping RADAR-IDs to MP-IDs.\nDefaults to [datadir]/[datadir].csv\n')
    cmdline.add_argument('--only-active', help='Only process subjects marked as active.\n', action="store_true")
    args = cmdline.parse_args()

    print("Processing directories {}".format(args.datadirs))

    for d in [ os.path.abspath(d) for d in args.datadirs]:
        study_id = os.path.basename(d)
        print()
        print("Study ID: {}".format(study_id))
        print("Calculating coverage report (all data)")

        if not os.path.isdir(d):
            print("{} is not a directory that exists, skipping.".format(d))
            continue

        # get patient mapping file path
        if not args.patient_map:
            map_path = os.path.join(d, study_id+".csv")
        else:
            map_path = args.patient_map

        if args.patient_ids:
            patients = [ {'RADAR-ID':pid, 'MP-ID':pid, 'status':'active'} for pid in args.patient_ids ]
        elif not os.path.isfile(map_path):
            print("No ID mapping file found for [{}], skipping.".format(d))
            continue
        else:
            # read patient mapping file
            patients = []
            with open(map_path) as map_csvfile:
                map_csvreader = csv.DictReader(map_csvfile)
                patients = [ p for p in map_csvreader ]

        coverage_all = {} # overall coverage store
        num_samples_all = {} # number of samples recorded
        for p in patients:
            if args.only_active and p["status"] != "active": continue
            coverage_all[p["MP-ID"]] = {} # initialize empty
            num_samples_all[p["MP-ID"]] = {} # initialize empty

            if not os.path.isdir(os.path.join(d, p["MP-ID"])):
                print("Data directory for {} ({}) not found, skipping.".format(p["RADAR-ID"], p["MP-ID"]))
                continue
            else:
                print("Processing {} ({})".format(p["RADAR-ID"], p["MP-ID"]))

            # get the first and last day where data is available
            (dt_first_day, dt_last_day) = getFirstLastDay({**COVERAGE_TOPICS, **QUESTIONNAIRE_TOPICS}, os.path.join(d, p["MP-ID"]))
            if any(i is None for i in [dt_first_day, dt_last_day]):
                print("No data available, skipping.")
                continue

            # calculate coverage value for each topic
            for t in COVERAGE_TOPICS:
                topic_coverage_all = [np.nan] # initialize empty
                topic_num_samples_all = 0 # initialize empty

                inpath = os.path.join(d, p["MP-ID"], t, "*.csv*")
                expected_lines = 3600 * COVERAGE_TOPICS[t] # 1h in sec * sample rate
                (topic_coverage_all, topic_num_samples_all) = getTopicCoverage(inpath, expected_lines, dt_first_day, dt_last_day)

                coverage_all[p["MP-ID"]][t] = topic_coverage_all
                num_samples_all[p["MP-ID"]][t] = topic_num_samples_all

        hours_rec = []
        for pid in num_samples_all.keys():
            hours_rec.append({"id":pid})
            for t in num_samples_all[pid]:
                hours = num_samples_all[pid][t] / COVERAGE_TOPICS[t] / 3600
                hours_rec[-1].update({getShortTopicStr(t):round(hours,2)})
        # pprint(hours_rec)

        days_rec = []
        for pid in coverage_all.keys():
            days_rec.append({"id":pid})
            t = 'android_empatica_e4_acceleration'

            days = set([ k[:8] for k in coverage_all[pid][t].keys() ])
            perday = { d:np.mean([ v for k,v in coverage_all[pid][t].items() if d in k ]) for d in days }
            #pprint(perday)
            dayslt5 = len([ v for v in perday.values() if v <= 0.05 ]) # days < 5%
            dayslt50 = len([ v for v in perday.values() if v > 0.05 and v <= 0.5 ]) # days < 50%
            daysgt50 = len([ v for v in perday.values() if v > 0.5 ]) # days > 50%
            days_rec[-1].update({'dayslt5':dayslt5,'dayslt50':dayslt50,'daysgt50':daysgt50})
        #print(days_rec)

        print('')
        print('Hours/Seizures Recorded:')
        fieldnames = ["id"] + [ getShortTopicStr(t) for t in COVERAGE_TOPICS.keys() ]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for row in hours_rec:
            writer.writerow(row)

        print('')
        print('Days Recorded:')
        fieldnames = [ 'id','dayslt5','dayslt50','daysgt50' ]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for row in days_rec:
            writer.writerow(row)
