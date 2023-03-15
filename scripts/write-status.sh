#!/bin/bash

current_date=$(date +"%Y-%m-%d %T") 
echo $current_date

if [ -z "$3" ];
	then 
		echo $current_date $1 >> $2;
	else
		FILE=$4
		if test -f "$FILE"; 
		then 
			echo $current_date $1 >> $2;
		else 
			echo "$current_date RaSP calculation failed (most likely the chain consisted only of hetatms)" >> $2;
		fi
	fi
