#!/usr/bin/env python3

import sys, os, glob, copy
import argparse, fileinput
import csv, gzip, json
from datetime import datetime, timedelta
import pytz
import math, re
import numpy as np
import matplotlib
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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



def plotCoveragePerday(coverage_data, figtitle, savepath=None):
    fig, ax = plt.subplots(figsize=(19,10))
    plt.set_cmap('RdYlGn')
    fig.suptitle(figtitle, fontsize='x-large')

    num_pat = len(coverage_data.keys())
    num_days = max([len(v) for v in coverage_data.values()])+180
    cov_mat = np.zeros((num_pat,num_days))
    for key_ix, key in enumerate(sorted(coverage_data.keys())):
        for val_ix, val in enumerate(coverage_data[key]):
            cov_mat[key_ix][val_ix] = val

    im = ax.imshow(cov_mat, vmin=0, vmax=1, aspect='auto', extent=[1,num_days+1,num_pat+1,1])

    xtick_range = np.array(range(0,num_days,30))
    xtick_range[0] += 1
    ax.set_xticks(xtick_range+0.5)
    ax.set_xticklabels(xtick_range)
    ax.set_xlabel("Day of recording")
    ax.set_yticks(np.array(range(num_pat))+1.5)
    ax.set_yticklabels(sorted(coverage_data.keys()))
    ax.set_ylim([num_pat+1,1])
    ax.set_ylabel("Participant ID")

    for d in range(2,num_days+1): ax.axvline(d, color='k', lw=1)
    for p in range(1,num_pat+2): ax.axhline(p, color='w', lw=5)
    for d in xtick_range[1:]: ax.axvline(d+1, color='k', lw=2)

    cb_ax = fig.add_axes([0.95, 0.2, 0.01, 0.6])
    cbar = fig.colorbar(im, cax=cb_ax, ticks=[0,0.25,0.5,0.75,1])

    if savepath:
        if not os.path.isdir(os.path.dirname(savepath)): os.mkdir(os.path.dirname(savepath))
        plt.savefig(savepath)
    else:
        plt.show()
    
    plt.close(fig)
    del fig
    plt.close('all')



if __name__=="__main__":
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass
    cmdline = argparse.ArgumentParser(description="RADAR-CNS data overview", formatter_class=Formatter)
    cmdline.add_argument('datadirs', type=str, nargs='+', help="data directories\n") # input directories
    cmdline.add_argument('--patient-ids', type=str, nargs='*', help='List of patient IDs (i.e. folder names), takes precedent over \'--patient-map\'.\n')
    cmdline.add_argument('--patient-map', metavar="PATH", type=str, help='Path to CSV file mapping RADAR-IDs to MP-IDs.\nDefaults to [datadir]/[datadir].csv\n')
    cmdline.add_argument('--only-active', help='Only process subjects marked as active.\n', action="store_true")
    cmdline.add_argument('-p', '--plot', help='Plot completeness per day per participant.\n', action="store_true")
    cmdline.add_argument('-o', '--outdir', type=str, default="plots", help="output directory\n")
    cmdline.add_argument('-f', '--outfmt', type=str, default="png", help="output format, i.e. file extension\n")
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

        days_rec = []
        for pid in coverage_all.keys():
            days_rec.append({"id":pid})
            t = 'android_empatica_e4_acceleration'

            days = set([ k[:8] for k in coverage_all[pid][t].keys() ])
            perday = { d:np.mean([ v for k,v in coverage_all[pid][t].items() if d in k ]) for d in days }
            dayslt5 = len([ v for v in perday.values() if v <= 0.05 ]) # days < 5%
            dayslt50 = len([ v for v in perday.values() if v > 0.05 and v <= 0.5 ]) # days < 50%
            daysgt50 = len([ v for v in perday.values() if v > 0.5 ]) # days > 50%
            days_rec[-1].update({'dayslt5':dayslt5,'dayslt50':dayslt50,'daysgt50':daysgt50,'perday':perday})

        if args.plot:
            # repack per day coverage data
            coverage_perday = { patcomp['id']:[patcomp['perday'][k] for k in sorted(patcomp['perday'].keys())] for patcomp in days_rec }

            # need to map from MP-ID to RADAR-ID here so it will show these in the plot correctly
            # first assert if RADAR-IDs are also unique, doesn't work otherwise
            radarIDs = [ pat["RADAR-ID"] for pat in patients ]
            assert len(radarIDs) == len(set(radarIDs))
            for oldkey in list(coverage_perday.keys()):
                newkey = next((item for item in patients if item["MP-ID"] == oldkey), {"RADAR-ID":"N/A"})["RADAR-ID"]
                coverage_perday[newkey] = coverage_perday.pop(oldkey)

            # call plotter
            figtitle = "Daily data coverage percentage"
            if args.outdir: savepath = os.path.abspath(os.path.join(args.outdir, "{}_dailyDataCoverage_{}.{}".format(study_id, datetime.now().strftime('%Y%m%d'), args.outfmt)))
            plotCoveragePerday(coverage_perday, figtitle, savepath if args.outdir else None)

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
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in days_rec:
            writer.writerow(row)
