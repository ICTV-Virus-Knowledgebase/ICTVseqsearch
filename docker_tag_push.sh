#!/usr/bin/env bash
#
# tag and push our latest image to docker hub
#
PUSH="no"
if [ "-push" == "$1" ]; then
    PUSH="yes"
fi
echo "PUSH=$PUSH"

# get version+git_hash
VER_FILE=version_git.txt
if [ ! -e "$VER_FILE" ]; then 
	echo "ERROR: $VER_FILE not found"
	exit 1
fi
VER_STR=$(cat $VER_FILE)


#
# tag
# 
echo sudo docker tag ictv_sequence_classifier curtish/ictv_sequence_classifier:$VER_STR
     sudo docker tag ictv_sequence_classifier curtish/ictv_sequence_classifier:$VER_STR
 
echo sudo docker tag ictv_sequence_classifier curtish/ictv_sequence_classifier:latest
     sudo docker tag ictv_sequence_classifier curtish/ictv_sequence_classifier:latest

#
# tag
#
if [ "$PUSH" == "yes" ]; then echo -n "# "; fi
echo sudo docker push curtish/ictv_sequence_classifier:$VER_STR
if [ "$PUSH" == "yes" ]; then
    sudo docker push curtish/ictv_sequence_classifier:$VER_STR
fi

if [ "$PUSH" == "yes" ]; then echo -n "# "; fi
echo sudo docker push curtish/ictv_sequence_classifier:latest
if [ "$PUSH" == "yes" ]; then
     sudo docker push curtish/ictv_sequence_classifier:latest
fi

#
# to pull instructions
#
echo "# to pull, run:"
echo "docker pull curtish/ictv_sequence_classifier:$VER_STR"
echo "docker pull curtish/ictv_sequence_classifier:latest"
