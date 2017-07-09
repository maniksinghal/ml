#!/usr/bin/python
import re
import sys
import subprocess
from os import listdir
from os.path import isfile, join
from datetime import datetime
import time
import argparse
import random

#Configuration
conf_ignore_duplicates = True

ENTRY_DELAY = 5000  # 5 seconds
DISCARD_ENTRY = "UNKNOWN"
g_features = list()
g_entry_db = list()

#
# Return absolute timestamp of the string (retrieved from ltrace)
# Or timediff between the string and previously calculated timestamp
# Returns value/diff in ms
def get_abs_time_diff_ms(line, old_time=None):
    #Jun 17 23:17:44.286 vic/evt 0/0/CPU0 t5361 0x0603 line 9.......
    match = re.search("(\S+\s+\d+\s+\d+:\d+:\d+)\.(\d+)\s+", line)
    if match != None:
	date_time = "1970 %s" % match.group(1)
	pattern = '%Y %b %d %H:%M:%S'
	value = datetime.strptime(date_time, pattern).timetuple()
        time_value = time.mktime(value)
	time_value = (time_value * 1000) + int(match.group(2))

	if old_time != None:
            return time_value - old_time
        else:
            return time_value
		
def process_file(file_path, port_filter):
    
    fd = open(file_path, "r")
    last_timestamp = 0
    event_ts = 0
    entry_features = None   # Features in current entry
    entry = None            # The actual entry (features + their value)
    link_status = DISCARD_ENTRY
    global g_entry_db 
    for line in fd:
        match = re.search("ether_ctrl_(\w+)_send_pulse.*%s.*fault=0x(\d).*\s(\w+)\svalue" % port_filter, line)
        if match == None:
            #Jun 17 18:29:48.076 vic/evt 0/0/CPU0 t5624 0x0410 line 574 vic_link_status_update( bay:0 if_inst:0 port:0 new_link_up=0x1 link_up=0x0 link_status_cnt=0x4d1 )
            match = re.search("vic_link_status_update.*%s.*new_link_up=0x(\d) link_up=0x(\d)" % port_filter, line)
            if match != None:
                if match.group(1) != match.group(2):
                    # link status changed
                    link_status = match.group(1)   # 1:UP, 0:Down
                        
            continue

        #print ("Match at: %s\n" % line)
        #Check timestamp first to decide whether its continuation of
        # current observation, or a new one
        ts = get_abs_time_diff_ms(line)
        ts = ts - last_timestamp
        last_timestamp = last_timestamp + ts
        if (ts > ENTRY_DELAY):

            # Store previous entry
            if entry != None and len(entry) != 0 and link_status != DISCARD_ENTRY:

                #entry.append(link_status)
                entry['link_status'] = link_status
                g_entry_db.append(entry)
                link_status = DISCARD_ENTRY

            entry_features = list()
            event_ts = 0.1  #First event marked as arrived at 0.1

            entry = {}
            entry_match = re.search("(\S+) vic/", line)
            if entry_match == None:
                print("Could not find string timestamp in %s" % line)
            entry['start_timestamp'] = entry_match.group(1)
            #entry.append(entry_match.group(1))  # Entry start timestamp in string format
            
        elif ts == 0:
            # Next event came at same time
            event_ts = event_ts + 0.1  # Mark next event as arrived 0.1ms later
        else:
            #Non-zero ts within the ENTRY_DELAY
            event_ts = event_ts + ts


        #Create feature variable
        pulse_type = match.group(1)
        fault = match.group(2)
        alarm = match.group(3)

        if fault == "1":
            fault = "ON"
        else:
            fault = "OFF"

        feature = "%s_%s_%s" % (alarm, pulse_type, fault)

        if conf_ignore_duplicates == False:
            i = 0
            while i < 20:
                if "%s_%d" % (feature, i) in entry_features:
                    i = i + 1
                else:
                    break
            feature = "%s_%d" % (feature, i)
            entry_features.append(feature)
            #entry.append("%s:%s" % (feature, event_ts))
            entry['%s' % feature] = event_ts
        else:
            if not feature in entry_features:
                entry_features.append(feature)
                #entry.append("%s:%s" % (feature, event_ts))
                entry['%s' % feature] = event_ts

        if not feature in g_features:
            g_features.append(feature)



    fd.close()
    return

def process_directory(path, port_filter):
    for f in listdir(path):
        file_path = join(path,f)
        if isfile(file_path):
            file_type = subprocess.check_output(['file', file_path])

            #Beware of tar/zip archives present in the directory
            if "ASCII text" in file_type:
                #print ("Processing %s" % file_path)
                process_file(file_path, port_filter)
                #print "REMOVE THIS RETURN\n"
                #return
            else:
                print("Ignored %s, as its not ASCII" % file_path)
    

def create_port_filter(interface, board):
    int_match = re.search("^\w+/(\d+)/(\d+)/(\d+)", interface)
    if int_match == None:
        print("Unknown interface type: %s" % interface)
        sys.exit()

    #@todo: Not considering board for now
    #@todo: Correct port calculation
    return "bay:%d.*port:%d" % (int(int_match.group(2)), int(int_match.group(3)))

parser = argparse.ArgumentParser(description="Model ltrace data for link flaps")
#parser.add_argument('--ws', default=None, help='Workspace to use')
parser.add_argument('interface', help="Interface to monitor")
parser.add_argument('board', help="Board type")
parser.add_argument('path', help="Directory path where ltraces are stored")

args = parser.parse_args()

#Create pattern to filter based on interface/board
bay_port_filter = create_port_filter(args.interface, args.board)
print("Created port-filter: %s" % bay_port_filter)

process_directory(args.path, bay_port_filter)
print ("Total features: %s\n" % g_features)


#Generate output
uniq = int(random.random() * 1000000)

fp = open('training_data_%d' % uniq, 'w')
mdata = open('training_Mdata_%d' % uniq, 'w')
#Print output
mdata.write("Label: 1(UP)/0(Down)\n")
for feature in g_features:
    mdata.write("%s\n" % feature)
mdata.write("\n\n")

for entry in g_entry_db:
    if entry['link_status'] == DISCARD_ENTRY:
        continue
    fp.write("%s" % entry['link_status'])
    for feature in g_features:
        if not feature in entry:
            fp.write(" 0")
        else:
            fp.write(" %s" % entry[feature])
    fp.write('\n')
fp.close()
mdata.close()

    

	
