#!/bin/bash

TGZ_FOLDER="$1"
TARGET_FOLDER="elk"

mkdir "$TARGET_FOLDER"
cd "$TARGET_FOLDER"
git init
cp $TGZ_FOLDER/README .
git add .
git commit -m "initial commit"

shopt -s nullglob
for f in ${TGZ_FOLDER}*.tgz
do
	# info
	echo "Adding file: $f"
	
	# clean current folder
	rm * -rf
	
	# unpack
	tar xzf $f --strip 1
	
	# git
	file_name=$(basename $f)
	git add .
	git commit -m "$file_name" --author="Elk <elk@elk.sourceforge.net>"
	git tag "$file_name"
done

