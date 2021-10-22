#!/bin/bash
# Compile data coverage report and send email.

cd "$( dirname "${BASH_SOURCE[0]}" )/.."

stack=bin/radar-docker
. lib/util.sh
. ./.env

#display_host="${SERVER_NAME} ($(hostname -f), $(curl -s http://ipecho.net/plain))"
display_host="${SERVER_NAME} ($(hostname -f), $(ip -f inet addr show ens5 | grep -Po 'inet \K[\d.]+'))"

function send_mail() {
    # Send mail notification to MAINTAINER
    # start up the mail container if not already started
    sudo-linux $stack up -d smtp

    # save the container, so that we can use exec to send an email later
    container=$(sudo-linux $stack ps -q smtp)
    address=$1
    subject=$2
    body=$3
    echo "Sending email notification to $address"
    echo -e "$body\nEmail sent: $(date -Iseconds)" | sudo-linux docker exec -i ${container} mail -aFrom:$FROM_EMAIL "-s${subject}" $address
}

echo "Checking incoming data and creating report."

startdt="$(date -Iseconds)"
report=$(sudo-linux python3 /home/ubuntu/radar-docker/lib/radar-data-coverage.py /mnt/hdfs-extract/UKLFR-Epi2 /mnt/hdfs-extract/KCL-Epi2)
echo -e "$report" | grep -e "-->" >> "/home/ubuntu/radar-docker/logs/radar-data-monitor/radar-data-coverage-$(date +%Gw%V).log"
enddt="$(date -Iseconds)"

subject="[RADAR-EPI-2] data coverage report for ${SERVER_NAME} ($(date -I))"
body="Data coverage report log for $(date):\n\n${report}\n\n\n----------\nScript start: ${startdt}"

#echo "$subject"
#echo "$body"

send_mail "mail@example.com" "$subject" "$body"

exit 0

