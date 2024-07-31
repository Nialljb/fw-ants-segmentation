#!/usr/bin/env bash 

# This script is used to run the gear locally for debugging purposes.
# It mounts the local directories to the docker container and runs the gear.
# The gear is run in the bash shell so that you can interact with the container.
# Assumes that API_KEY is set in the environment and added to config.json.

GEAR=fw-ants-segmentation
IMAGE=flywheel/ants-segmentation:0.2.4
LOG=ants-segmentation-0.2.4-66a7e400e1bb69b99fc4ae86

# Command:
docker run -it --cpus 6.0 --rm --entrypoint bash\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/app/:/flywheel/v0/app\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/utils:/flywheel/v0/utils\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/run.py:/flywheel/v0/run.py\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/${LOG}/input:/flywheel/v0/input\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/${LOG}/output:/flywheel/v0/output\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/${LOG}/work:/flywheel/v0/work\
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/${LOG}/config.json:/flywheel/v0/config.json \
	-v /Users/nbourke/GD/atom/unity/fw-gears/${GEAR}/${LOG}/manifest.json:/flywheel/v0/manifest.json\
	$IMAGE
