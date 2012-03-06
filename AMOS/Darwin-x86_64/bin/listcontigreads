#!/bin/sh

grep '#' $1 | grep -v '##' | cut -f1 -d')' | tr -d '#' | tr '(' ' '  | sort -nk2 | awk '{print $1}'
