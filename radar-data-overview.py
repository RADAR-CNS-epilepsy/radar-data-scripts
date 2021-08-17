#!/usr/bin/env python3

import sys, os, glob, copy
import argparse, fileinput
import csv, gzip, json
from datetime import datetime, timedelta
import pytz
import math, re
import numpy as np
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
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

def plotCoverage(coverage_data, figtitle, savepath=None, q_dates=None, mpid = "N/A"):
    fig, axs = plt.subplots(len(coverage_data), 1, figsize=(19, 10))
    plt.set_cmap('RdYlGn')
    fig.suptitle(figtitle, fontsize='x-large')

    for t_ix, t in enumerate(coverage_data.keys()):
        days = sorted(set([ datetime.strptime(d, FILEMASK).replace(hour=0) for d in sorted(coverage_data[t].keys()) ]))
        num_days = int(np.ceil(len(coverage_data[t])/24))
        assert len(days) == num_days
        cov_mat = np.empty((24,num_days))
        cov_mat.fill(np.nan)

        for key_ix, h in enumerate(sorted(coverage_data[t].keys())):
            h_ix = datetime.strptime(h, FILEMASK).hour
            d_ix = int(np.floor(key_ix/24))
            cov_mat[h_ix,d_ix] = coverage_data[t][h]

        axs[t_ix].set_title("-".join(t.split('_')[1:]))
        axs[t_ix].set_ylim([24,0])
        im = axs[t_ix].imshow(cov_mat, aspect='auto', extent=[mdates.date2num(days[0]-timedelta(hours=12)), mdates.date2num(days[-1]+timedelta(hours=12)), 24, 0], vmin=0, vmax=1)
        axs[t_ix].xaxis_date()
        axs[t_ix].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        axs[t_ix].xaxis.set_major_locator(mdates.DayLocator(interval=int(np.ceil(num_days/30))))
        for d in range(num_days-1): axs[t_ix].axvline(mdates.date2num(days[0] + timedelta(days=d, hours=12)), color='k')
        axs[t_ix].set_xlabel("day of recording")
        axs[t_ix].set_ylabel("hour of day")

    fig.autofmt_xdate()
    # Make space for title
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.05, right=0.95, hspace=0.25)
    cb_ax = fig.add_axes([0.96, 0.2, 0.01, 0.6])
    cbar = fig.colorbar(im, cax=cb_ax)

    # mark questionnaire dates
    if q_dates:
        # iterate topics
        markers, marker_labels = [],[] # for legend
        for q_topic_ix, q_topic in enumerate(q_dates):
            # calculate marker coordinates: date as mdate for x-axis, time of day as float for y-axis; only mark if within already plotted dates
            q_coords = np.array([ (mdates.date2num(q.date()), q.hour+q.minute/60.0) for q in q_dates[q_topic] if q.date() in [d.date() for d in days] ])
            # test if any to be marked
            if len(q_coords.shape) != 2: continue
            # mark in all subplots
            for ax in axs: marker = ax.scatter(q_coords[:,0], q_coords[:,1], s=(50 if num_days < 90 else 25), c=QUESTIONNAIRE_TOPICS[q_topic][1], marker=QUESTIONNAIRE_TOPICS[q_topic][0])
            # save info for legend
            markers.append(marker)
            marker_labels.append("-".join(q_topic.split('_')[1:]) + "; N={}".format(q_coords.shape[0]))

        fig.legend(markers, marker_labels, facecolor='lightgrey')

    fig.text(0.01, 0.98, "created: {}".format(datetime.now(pytz.utc).strftime("%F %T %z")))
    fig.text(0.01, 0.96, "MP-ID: {}".format(mpid))

    if savepath:
        if not os.path.isdir(os.path.dirname(savepath)): os.mkdir(os.path.dirname(savepath))
        plt.savefig(savepath)
    else:
        plt.show()
    
    plt.close(fig)
    del fig
    plt.close('all')


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


def getTopicCoverage(inpath, expected_lines, dt_first_day, dt_last_day, lookback_days):
    # init, return if no files found
    topic_coverage_all = [np.nan]
    topic_coverage_lookback = [np.nan]
    topic_num_samples_all = 0
    infiles = sorted(glob.glob(inpath))
    if not infiles: return (topic_coverage_all, topic_coverage_lookback, topic_num_samples_all)

    # create the coverage value timelines
    dt_today = datetime.utcnow().replace(minute=0, second=0, microsecond=0, hour=23)
    num_hours_all = int((dt_last_day-dt_first_day).total_seconds()/3600)
    expected_hours_all = [ dt_first_day + x * HOUR for x in range(num_hours_all+1) ]
    expected_hours_lookback = [ dt_today - x * HOUR for x in range(lookback_days*24) ]
    topic_coverage_all = {datetime.strftime(h, FILEMASK):0 for h in expected_hours_all}
    topic_coverage_lookback = {datetime.strftime(h, FILEMASK):0 for h in expected_hours_lookback}

    # iterate data files per topic, fill timeline with data coverage from available files
    for f in infiles:
        fsplit = os.path.basename(f).split('.')
        num_lines = 0

        # read lines, skip if not in timeline
        if fsplit[0] not in topic_coverage_all.keys() and fsplit[0] not in topic_coverage_lookback.keys(): continue
        if fsplit[-1] == "gz": num_lines = sum(1 for line in gzip.open(f))
        elif fsplit[-1] == "csv": num_lines = sum(1 for line in open(f))
        else: continue

        # calculate coverage, fill in timeline, upper bound at 100%
        file_coverage = num_lines/expected_lines
        if fsplit[0] in topic_coverage_all.keys(): topic_coverage_all[fsplit[0]] = file_coverage if file_coverage < 1 else 1
        if fsplit[0] in topic_coverage_lookback.keys(): topic_coverage_lookback[fsplit[0]] = file_coverage if file_coverage < 1 else 1
        topic_num_samples_all += num_lines
    
    return (topic_coverage_all, topic_coverage_lookback, topic_num_samples_all)


def getQuestionnaireCompletedDates(topics, inpath, datefield):
    q_dates = {}
    for t in topics:
        q_dates[t] = []
        # get all the csv files in the topic folder
        infiles = sorted(glob.glob(os.path.join(inpath, t, "*.csv*")))
        if not infiles: continue

        # read questionnaire csv files
        q_events = []
        for q_infile in infiles:
            q_fsplit = os.path.basename(q_infile).split('.')
            if q_fsplit[-1] == "gz":
                with gzip.open(q_infile, 'rt') as q_csvfile:
                    q_csvreader = csv.DictReader(q_csvfile)
                    q_events.extend([ q for q in q_csvreader ])
            elif q_fsplit[-1] == "csv":
                with open(q_infile) as q_csvfile:
                    q_csvreader = csv.DictReader(q_csvfile)
                    q_events.extend([ q for q in q_csvreader ])
            else: continue

        # get completed times
        for qe in q_events:
            if "seizure_diary" in t:
                sz_stamp = getUTCFromSDAnswer(qe[datefield], guessTZFromPath(inpath))
                q_dates[t].append(datetime.fromtimestamp(sz_stamp,tz=pytz.utc))
            else:
                q_dates[t].append(datetime.fromtimestamp(float(qe[datefield]),tz=pytz.utc))
    return q_dates


def getMonthFromStr(mstr):
    if mstr == 'Jan' or mstr == 'Jan.': return 1
    elif mstr == 'Feb' or mstr == 'Feb.': return 2
    elif mstr == 'Mar' or mstr == 'März': return 3
    elif mstr == 'Apr' or mstr == 'Apr.': return 4
    elif mstr == 'May' or mstr == 'Mai': return 5
    elif mstr == 'Jun' or mstr == 'Juni': return 6
    elif mstr == 'Jul' or mstr == 'Juli': return 7
    elif mstr == 'Aug' or mstr == 'Aug.': return 8
    elif mstr == 'Sep' or mstr == 'Sep.': return 9
    elif mstr == 'Oct' or mstr == 'Okt.': return 10
    elif mstr == 'Nov' or mstr == 'Nov.': return 11
    elif mstr == 'Dec' or mstr == 'Dez.': return 12
    else: return -1
def getHoursFromAMPM(hstr,ampm):
    h = int(hstr)
    if h == 12 and 'AM' in ampm: h = 0
    elif h != 12 and 'PM' in ampm: h = h + 12
    return h
def guessTZFromPath(path):
    if "UKF" in path or "UKLFR" in path: return "Europe/Berlin"
    elif "KCL" in path: return "Europe/London"
    else: return "Etc/UCT"
def getUTCFromSDAnswer(sd_value, tz):
    sd_info = json.loads(sd_value)
    dt_info = datetime(year=int(sd_info['year']), month=getMonthFromStr(sd_info['month']), day=int(sd_info['day']), hour=getHoursFromAMPM(sd_info['hour'],sd_info['ampm']), minute=int(sd_info['minute']), second=int(sd_info['second']))
    dt_info = pytz.timezone(tz).localize(dt_info)
    return dt_info.timestamp()


def getNumSeizuresRecorded(sz_list, inpath, margin_min=10):
    infiles = sorted(glob.glob(inpath))
    if not infiles: return 0

    num_sz_recorded = 0
    for sz in sz_list:
        # get relevant data files
        sz_range = [sz - timedelta(minutes=margin_min), sz , sz + timedelta(minutes=margin_min)]
        sz_filenames = [ dt.strftime('%Y%m%d_%H00.csv') for dt in sz_range ]
        matches = [ f for f in infiles if any([re.search(sz_f, f) is not None for sz_f in sz_filenames]) ]

        # read data
        data = []
        for m in matches:
            fsplit = os.path.basename(m).split('.')
            if fsplit[-1] == "gz": openfn = gzip.open
            elif fsplit[-1] == "csv": openfn = open
            else: continue
            with openfn(m, mode='rt') as csvfile:
                csvreader = csv.DictReader(csvfile)
                for row in csvreader: data.append(row)

        # assess seizure data completeness
        sz_min = sz_range[0].timestamp()
        sz_max = sz_range[2].timestamp()
        data_samples_dt = np.array([int(float(d['value.time'])) for d in data])
        sz_samples_dt = np.where((data_samples_dt > sz_min) & (data_samples_dt < sz_max))[0]
        sz_data_th = 0.9*(2*margin_min)*(60*32) # 90% of 2*margin_min minutes ACC data
        if sz_samples_dt.size >= sz_data_th: num_sz_recorded += 1

    return num_sz_recorded




if __name__=="__main__":
    class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter): pass
    cmdline = argparse.ArgumentParser(description="RADAR-CNS data overview", formatter_class=Formatter)
    cmdline.add_argument('datadirs', type=str, nargs='+', help="data directories\n") # input directories
    cmdline.add_argument('--patient-ids', type=str, nargs='*', help='List of patient IDs (i.e. folder names), takes precedent over \'--patient-map\'.\n')
    cmdline.add_argument('--patient-map', metavar="PATH", type=str, help='Path to CSV file mapping RADAR-IDs to MP-IDs.\nDefaults to [datadir]/[datadir].csv\n')
    cmdline.add_argument('-l', '--lookback', metavar="DAYS", type=int, default=14, help='Number of days to look back for coverage report.\n')
    cmdline.add_argument('--only-active', help='Only process subjects marked as active.\n', action="store_true")
    cmdline.add_argument('-o', '--outdir', type=str, default="plots", help="output directory\n")
    cmdline.add_argument('-f', '--outfmt', type=str, default="png", help="output format, i.e. file extension\n")
    cmdline.add_argument('-np', '--no-plot', help='Do not plot, only print log to stdout.\n', action="store_true")
    args = cmdline.parse_args()

    print("Processing directories {}".format(args.datadirs))

    for d in [ os.path.abspath(d) for d in args.datadirs]:
        study_id = os.path.basename(d)
        print()
        print("Study ID: {}".format(study_id))
        print("Plotting data coverage report (all data + last {} days)".format(args.lookback))

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
        coverage_lookback = {} # lookback coverage store
        num_samples_all = {} # number of samples recorded
        num_seizures_all = {} # number of seizures reported/recorded
        for p in patients:
            if args.only_active and p["status"] != "active": continue
            coverage_all[p["MP-ID"]] = {} # initialize empty
            coverage_lookback[p["MP-ID"]] = {} # initialize empty
            num_samples_all[p["MP-ID"]] = {} # initialize empty
            num_seizures_all[p["MP-ID"]] = {} # initialize empty

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
                topic_coverage_lookback = [np.nan] # initialize empty
                topic_num_samples_all = 0 # initialize empty

                inpath = os.path.join(d, p["MP-ID"], t, "*.csv*")
                expected_lines = 3600 * COVERAGE_TOPICS[t] # 1h in sec * sample rate
                (topic_coverage_all, topic_coverage_lookback, topic_num_samples_all) = getTopicCoverage(inpath, expected_lines, dt_first_day, dt_last_day, args.lookback)

                coverage_all[p["MP-ID"]][t] = topic_coverage_all
                coverage_lookback[p["MP-ID"]][t] = topic_coverage_lookback
                num_samples_all[p["MP-ID"]][t] = topic_num_samples_all

            # get questionnaire completed dates
            #q_dates = getQuestionnaireCompletedDates(QUESTIONNAIRE_TOPICS, os.path.join(d, p["MP-ID"]), "value.timeCompleted")
            q_dates = getQuestionnaireCompletedDates(['questionnaire_epi_evening_questionnaire'], os.path.join(d, p["MP-ID"]), "value.timeCompleted")
            s_dates = getQuestionnaireCompletedDates(['questionnaire_epi_seizure_diary'], os.path.join(d, p["MP-ID"]), "value.answers.0.value")
            q_dates.update(s_dates)

            # get seizure reported/recorded numbers
            if 'questionnaire_epi_seizure_diary' in q_dates.keys():
                sz = q_dates['questionnaire_epi_seizure_diary']
                num_sz_rep = len(sz)
                num_sz_rec = getNumSeizuresRecorded(sz, os.path.join(d, p["MP-ID"], 'android_empatica_e4_acceleration', "*.csv*"))
                num_seizures_all[p["MP-ID"]] = {'sz_rep':num_sz_rep, 'sz_rec':num_sz_rec}

            # plotting
            if not args.no_plot:
                num_days = max([ int(np.ceil(len(coverage_all[p["MP-ID"]][t])/24)) for t in coverage_all[p["MP-ID"]].keys() ])
                figtitle_all = "hourly data coverage percentage for '{}' (all data, {} days)".format(p["RADAR-ID"], num_days)
                figtitle_lookback = "hourly data coverage percentage for '{}' (last {} days)".format(p["RADAR-ID"], args.lookback)
                if args.outdir:
                    savepath_all = os.path.abspath(os.path.join(args.outdir, "{}_{}_{}_covall.{}".format(study_id, p["RADAR-ID"], p["MP-ID"], args.outfmt)))
                    savepath_lookback = os.path.abspath(os.path.join(args.outdir, "{}_{}_{}_cov{}d.{}".format(study_id, p["RADAR-ID"], p["MP-ID"], args.lookback, args.outfmt)))
                plotCoverage(coverage_all[p["MP-ID"]], figtitle_all, savepath_all if args.outdir else None, q_dates, mpid=p["MP-ID"])
                plotCoverage(coverage_lookback[p["MP-ID"]], figtitle_lookback, savepath_lookback if args.outdir else None, q_dates, mpid=p["MP-ID"])

        hours_rec = []
        for pid in num_samples_all.keys():
            hours_rec.append({"id":pid})
            for t in num_samples_all[pid]:
                hours = num_samples_all[pid][t] / COVERAGE_TOPICS[t] / 3600
                hours_rec[-1].update({getShortTopicStr(t):round(hours,2)})
            hours_rec[-1].update(num_seizures_all[pid])
        #pprint(hours_rec)

        print('')
        print('Hours/Seizures Recorded:')
        fieldnames = ["id"] + [ getShortTopicStr(t) for t in COVERAGE_TOPICS.keys() ] + ["sz_rep", "sz_rec"]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        for row in hours_rec:
            writer.writerow(row)
