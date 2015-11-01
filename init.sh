#! /bin/bash

# Check for config file
if [ ! -e config.json ]; then
	printf "BASH: Fatal Error - Workflow file missing (config.json)"
	exit
fi

# Check for core python module
if [ ! -e core.py ]; then
	printf "BASH: Fatal Error - Workflow file missing (core.py)"
	exit
fi

printf "\nBASH: Status - File check successful, gathering permissions to execute workflow\n"

# Assign permissions to directory
{
	printf "\nBASH: Notification - Attempting recursive modification of assignments directory"
	printf "\n\t 'chown -R $(whoami):admin assignments'"
	printf "\n\t 'chmod -R 775 assignments'"
	chown -R $(whoami):admin assignments
	chmod -R 775 assignments
} || {
	printf "\nBASH: Notification - Unable to perform permissions modification; running as superuser - please enter your OS password\n"
	sudo chown -R $(whoami):admin assignments
	sudo chmod -R 775 assignments
}

printf "\nBASH: Status - Permissions modified successfully\n"
printf "\nBASH: Status - Executing core.py\n"

# Init python core module
python core.py