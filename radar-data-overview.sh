#!/bin/bash
# Plot data coverage overviews, aggregate to pdf, and send by email.

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
    sentstr="Email sent: $(date -Iseconds)"
    sendbody=$(echo -e "${body}\n${sentstr}")
    attachment=$4
    echo "Sending email notification to $address"
    mime-construct --output --to "$address" --subject "$subject" --encoding quoted-printable --string "${sendbody}" --file-attach ${attachment} | sudo-linux docker exec -i ${container} mail -t
}

echo "Processing extracted data and plotting coverage."

startdt="$(date -Iseconds)"
pdfout="/home/ubuntu/radar-docker/coverage-plots/coverage-plots-$(date -I).pdf"
report=$(sudo-linux python3 /home/ubuntu/radar-docker/lib/radar-data-overview.py /mnt/hdfs-extract/UKLFR-Epi2 /mnt/hdfs-extract/KCL-Epi2 --only-active -o /home/ubuntu/radar-docker/coverage-plots)
img2pdf /home/ubuntu/radar-docker/coverage-plots/*.png -o "$pdfout"
enddt="$(date -Iseconds)"

subject="[RADAR-EPI-2] data coverage weekly digest for ${SERVER_NAME} ($(date -I))"
body="Weekly digest of data coverage for $(date), see attachment.\n\n\n----------\nScript log:\n${report}\n\n----------\nScript start: ${startdt}"
attachment=${pdfout}

#echo "$subject"
#echo "$body"
#echo "$attachment"

send_mail "mail@example.com" "$subject" "$body" "$attachment"

exit 0

