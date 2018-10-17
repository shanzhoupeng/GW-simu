#!/bin/bash
code=$(curl 'http://localhost:17304/' | tr -cd "[0-9]")
code4=${code#*33}
./loginhpc $code4
