#!/bin/sh

submodule_init ( ) {
    local SUBMODULES=$(git submodule | awk '{print $2}')
    for submodule in $SUBMODULES; do
	echo "Initializing submodule $submodule"
	git submodule init $submodule
	git submodule update $submodule
	cd $submodule
	if test -f version.sh
	then
	    ./version.sh
	fi
	submodule_init
	cd ..
    done
}

./version.sh
submodule_init
autoreconf --install
