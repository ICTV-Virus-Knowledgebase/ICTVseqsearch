#!/usr/bin/env bash
#
#
if [ -z "$(which rclone 2>/dev/null)" ]; then
    echo module load rclone
    module load rclone
fi

SRC_URL="box:/Virus Knowledgebase/VMR-blast/VMRs"
DEST_URL="./VMRs/"

echo rclone --verbose copy "$SRC_URL" "$DEST_URL"
     rclone --verbose copy "$SRC_URL" "$DEST_URL"

