#!/bin/sh
# Reproduces every number and figure of the study (about 20 s).
set -e
cd "$(dirname "$0")/src"
python3 trigger_eq.py
python3 frontend250.py
python3 system.py
