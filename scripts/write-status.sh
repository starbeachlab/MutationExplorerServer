#!/bin/bash

current_date=$(date +"%Y-%m-%d %T") 
echo $current_date
echo $current_date $1 >> $2
